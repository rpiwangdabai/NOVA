/*
 * BRIEF DESCRIPTION
 *
 * DAX file operations.
 *
 * Copyright 2015-2016 Regents of the University of California,
 * UCSD Non-Volatile Systems Lab, Andiry Xu <jix024@cs.ucsd.edu>
 * Copyright 2012-2013 Intel Corporation
 * Copyright 2009-2011 Marco Stornelli <marco.stornelli@gmail.com>
 *
 * This file is licensed under the terms of the GNU General Public
 * License version 2. This program is licensed "as is" without any
 * warranty of any kind, whether express or implied.
 */

#include <linux/buffer_head.h>
#include <asm/cpufeature.h>
#include <asm/pgtable.h>
#include <linux/version.h>
#include "nova.h"
#include "journal.h"
#ifdef DYNAMIC_JOURNAL
static size_t write_len_average =1024;
static u64 write_count = 0;
static u64 write_len_sum = 0;
#endif
	static ssize_t
do_dax_mapping_read(struct file *filp, char __user *buf,
	size_t len, loff_t *ppos)
{
	struct inode *inode = filp->f_mapping->host;
	struct super_block *sb = inode->i_sb;
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	struct nova_file_write_entry *entry;
	pgoff_t index, end_index;
	unsigned long offset;
	loff_t isize, pos;
	size_t copied = 0, error = 0;
	timing_t memcpy_time;

	pos = *ppos;
	index = pos >> PAGE_SHIFT;
	offset = pos & ~PAGE_MASK;

	if (!access_ok(VERIFY_WRITE, buf, len)) {
		error = -EFAULT;
		goto out;
	}

	isize = i_size_read(inode);
	if (!isize)
		goto out;

	nova_dbgv("%s: inode %lu, offset %lld, count %lu, size %lld\n",
		__func__, inode->i_ino,	pos, len, isize);

	if (len > isize - pos)
		len = isize - pos;

	if (len <= 0)
		goto out;

	end_index = (isize - 1) >> PAGE_SHIFT;
	do {
		unsigned long nr, left;
		unsigned long nvmm;
		void *dax_mem = NULL;
		int zero = 0;

		/* nr is the maximum number of bytes to copy from this page */
		if (index >= end_index) {
			if (index > end_index)
				goto out;
			nr = ((isize - 1) & ~PAGE_MASK) + 1;
			if (nr <= offset) {
				goto out;
			}
		}

		entry = nova_get_write_entry(sb, si, index);
		if (unlikely(entry == NULL)) {
			nova_dbgv("Required extent not found: pgoff %lu, "
				"inode size %lld\n", index, isize);
			nr = PAGE_SIZE;
			zero = 1;
			goto memcpy;
		}

		/* Find contiguous blocks */
		if (index < entry->pgoff ||
			index - entry->pgoff >= entry->num_pages) {
			nova_err(sb, "%s ERROR: %lu, entry pgoff %llu, num %u, "
				"blocknr %llu\n", __func__, index, entry->pgoff,
				entry->num_pages, entry->block >> PAGE_SHIFT);
			return -EINVAL;
		}
		if (entry->invalid_pages == 0) {
			nr = (entry->num_pages - (index - entry->pgoff))
				* PAGE_SIZE;
		} else {
			nr = PAGE_SIZE;
		}

		nvmm = get_nvmm(sb, sih, entry, index);
		dax_mem = nova_get_block(sb, (nvmm << PAGE_SHIFT));

memcpy:
		nr = nr - offset;
		if (nr > len - copied)
			nr = len - copied;

		NOVA_START_TIMING(memcpy_r_nvmm_t, memcpy_time);

		if (!zero)
			left = __copy_to_user(buf + copied,
						dax_mem + offset, nr);
		else
			left = __clear_user(buf + copied, nr);

		NOVA_END_TIMING(memcpy_r_nvmm_t, memcpy_time);

		if (left) {
			nova_dbg("%s ERROR!: bytes %lu, left %lu\n",
				__func__, nr, left);
			error = -EFAULT;
			goto out;
		}

		copied += (nr - left);
		offset += (nr - left);
		index += offset >> PAGE_SHIFT;
		offset &= ~PAGE_MASK;
	} while (copied < len);

out:
	*ppos = pos + copied;
	if (filp)
		file_accessed(filp);

	NOVA_STATS_ADD(read_bytes, copied);

	nova_dbgv("%s returned %zu\n", __func__, copied);
	return (copied ? copied : error);
}

/*
 * Wrappers. We need to use the rcu read lock to avoid
 * concurrent truncate operation. No problem for write because we held
 * i_mutex.
 */
ssize_t nova_dax_file_read(struct file *filp, char __user *buf,
			    size_t len, loff_t *ppos)
{
	ssize_t res;
	timing_t dax_read_time;

	NOVA_START_TIMING(dax_read_t, dax_read_time);
//	rcu_read_lock();
	res = do_dax_mapping_read(filp, buf, len, ppos);
//	rcu_read_unlock();
	NOVA_END_TIMING(dax_read_t, dax_read_time);
	return res;
}

static inline int nova_copy_partial_block(struct super_block *sb,
	struct nova_inode_info_header *sih,
	struct nova_file_write_entry *entry, unsigned long index,
	size_t offset, void* kmem, bool is_end_blk)
{
	void *ptr;
	unsigned long nvmm;

	nvmm = get_nvmm(sb, sih, entry, index);
	ptr = nova_get_block(sb, (nvmm << PAGE_SHIFT));
	if (ptr != NULL) {
		if (is_end_blk) {
			memcpy(kmem + offset, ptr + offset,
				sb->s_blocksize - offset);
#ifdef NVM_DELAY
			ndelay((((sb->s_blocksize - offset-1)>>NVM_BLOCK_SHIFT)+1)<<BLOCK_DELAY_SHIFT);
#endif
		}
		else 
			memcpy(kmem, ptr, offset);
#ifdef NVM_DELAY
			ndelay((((offset-1)>>NVM_BLOCK_SHIFT)+1)<<BLOCK_DELAY_SHIFT);
#endif
	}

	return 0;
}

/* 
 * Fill the new start/end block from original blocks.
 * Do nothing if fully covered; copy if original blocks present;
 * Fill zero otherwise.
 */
static void nova_handle_head_tail_blocks(struct super_block *sb,
	struct nova_inode *pi, struct inode *inode, loff_t pos, size_t count,
	void *kmem)
{
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	size_t offset, eblk_offset;
	unsigned long start_blk, end_blk, num_blocks;
	struct nova_file_write_entry *entry;
	timing_t partial_time;

	NOVA_START_TIMING(partial_block_t, partial_time);
	offset = pos & (sb->s_blocksize - 1);
	num_blocks = ((count + offset - 1) >> sb->s_blocksize_bits) + 1;
	/* offset in the actual block size block */
	offset = pos & (nova_inode_blk_size(pi) - 1);
	start_blk = pos >> sb->s_blocksize_bits;
	end_blk = start_blk + num_blocks - 1;

	nova_dbg_verbose("%s: %lu blocks\n", __func__, num_blocks);
	/* We avoid zeroing the alloc'd range, which is going to be overwritten
	 * by this system call anyway */
	nova_dbg_verbose("%s: start offset %lu start blk %lu %p\n", __func__,
				offset, start_blk, kmem);
	if (offset != 0) {
		entry = nova_get_write_entry(sb, si, start_blk);
		if (entry == NULL) {
			/* Fill zero */
		    	memset(kmem, 0, offset);
		} else {
			/* Copy from original block */
			nova_copy_partial_block(sb, sih, entry, start_blk,
					offset, kmem, false);
		}
		nova_flush_buffer(kmem, offset, 0);
	}

	kmem = (void *)((char *)kmem +
			((num_blocks - 1) << sb->s_blocksize_bits));
	eblk_offset = (pos + count) & (nova_inode_blk_size(pi) - 1);
	nova_dbg_verbose("%s: end offset %lu, end blk %lu %p\n", __func__,
				eblk_offset, end_blk, kmem);
	if (eblk_offset != 0) {
		entry = nova_get_write_entry(sb, si, end_blk);
		if (entry == NULL) {
			/* Fill zero */
		    	memset(kmem + eblk_offset, 0,
					sb->s_blocksize - eblk_offset);
		} else {
			/* Copy from original block */
			nova_copy_partial_block(sb, sih, entry, end_blk,
					eblk_offset, kmem, true);
		}
		nova_flush_buffer(kmem + eblk_offset,
					sb->s_blocksize - eblk_offset, 0);
	}

	NOVA_END_TIMING(partial_block_t, partial_time);
}

int nova_reassign_file_tree(struct super_block *sb,
	struct nova_inode *pi, struct nova_inode_info_header *sih,
	u64 begin_tail)
{
	struct nova_file_write_entry *entry_data;
	u64 curr_p = begin_tail;
	size_t entry_size = sizeof(struct nova_file_write_entry);

	while (curr_p != pi->log_tail) {
		if (is_last_entry(curr_p, entry_size))
			curr_p = next_log_page(sb, curr_p);

		if (curr_p == 0) {
			nova_err(sb, "%s: File inode %llu log is NULL!\n",
				__func__, pi->nova_ino);
			return -EINVAL;
		}

		entry_data = (struct nova_file_write_entry *)
					nova_get_block(sb, curr_p);

		if (nova_get_entry_type(entry_data) != FILE_WRITE) {
			nova_dbg("%s: entry type is not write? %d\n",
				__func__, nova_get_entry_type(entry_data));
			curr_p += entry_size;
			continue;
		}

		nova_assign_write_entry(sb, pi, sih, entry_data, true);
		curr_p += entry_size;
	}

	return 0;
}


static int nova_cleanup_incomplete_write(struct super_block *sb,
	struct nova_inode *pi, struct nova_inode_info_header *sih,
	unsigned long blocknr, int allocated, u64 begin_tail, u64 end_tail)
{
	struct nova_file_write_entry *entry;
	u64 curr_p = begin_tail;
	size_t entry_size = sizeof(struct nova_file_write_entry);

	if (blocknr > 0 && allocated > 0)
		nova_free_data_blocks(sb, pi, blocknr, allocated);

	if (begin_tail == 0 || end_tail == 0)
		return 0;

	while (curr_p != end_tail) {
		if (is_last_entry(curr_p, entry_size))
			curr_p = next_log_page(sb, curr_p);

		if (curr_p == 0) {
			nova_err(sb, "%s: File inode %llu log is NULL!\n",
				__func__, pi->nova_ino);
			return -EINVAL;
		}

		entry = (struct nova_file_write_entry *)
					nova_get_block(sb, curr_p);

		if (nova_get_entry_type(entry) != FILE_WRITE) {
			nova_dbg("%s: entry type is not write? %d\n",
				__func__, nova_get_entry_type(entry));
			curr_p += entry_size;
			continue;
		}

		blocknr = entry->block >> PAGE_SHIFT;
		nova_free_data_blocks(sb, pi, blocknr, entry->num_pages);
		curr_p += entry_size;
	}

	return 0;
}


ssize_t nova_cow_file_write(struct file *filp,
	const char __user *buf,	size_t len, loff_t *ppos, bool need_mutex)
{
	struct address_space *mapping = filp->f_mapping;
	struct inode    *inode = mapping->host;
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	struct super_block *sb = inode->i_sb;
	struct nova_inode *pi;
	struct nova_file_write_entry entry_data;
	ssize_t     written = 0;
	loff_t pos;
	size_t count, offset, copied, ret;
	unsigned long start_blk, num_blocks;
	unsigned long total_blocks;
	unsigned long blocknr = 0;
	unsigned int data_bits;
	int allocated = 0;
	void* kmem;
	u64 curr_entry;
	size_t bytes;
	long status = 0;
	timing_t cow_write_time, memcpy_time;
	unsigned long step = 0;
	u64 temp_tail = 0, begin_tail = 0;
	u32 time;

	if (len == 0)
		return 0;

	/*
	 * We disallow writing to a mmaped file,
	 * since write is copy-on-write while mmap is DAX (in-place).
	 */
	if (mapping_mapped(mapping))
		return -EACCES;

	NOVA_START_TIMING(cow_write_t, cow_write_time);

	sb_start_write(inode->i_sb);
	if (need_mutex)
		mutex_lock(&inode->i_mutex);

	if (!access_ok(VERIFY_READ, buf, len)) {
		ret = -EFAULT;
		goto out;
	}
	pos = *ppos;

	if (filp->f_flags & O_APPEND)
		pos = i_size_read(inode);

	count = len;

	pi = nova_get_inode(sb, inode);

	offset = pos & (sb->s_blocksize - 1);
	num_blocks = ((count + offset - 1) >> sb->s_blocksize_bits) + 1;
	total_blocks = num_blocks;
	/* offset in the actual block size block */

	ret = file_remove_privs(filp);
	if (ret) {
		goto out;
	}
	inode->i_ctime = inode->i_mtime = CURRENT_TIME_SEC;
	time = CURRENT_TIME_SEC.tv_sec;

	nova_dbgv("%s: inode %lu, offset %lld, count %lu\n",
			__func__, inode->i_ino,	pos, count);

	temp_tail = pi->log_tail;
	while (num_blocks > 0) {
		offset = pos & (nova_inode_blk_size(pi) - 1);
		start_blk = pos >> sb->s_blocksize_bits;

		/* don't zero-out the allocated blocks */
		allocated = nova_new_data_blocks(sb, pi, &blocknr, num_blocks,
						start_blk, 0, 1);
		nova_dbg_verbose("%s: alloc %d blocks @ %lu\n", __func__,
						allocated, blocknr);

		if (allocated <= 0) {
			nova_dbg("%s alloc blocks failed %d\n", __func__,
								allocated);
			ret = allocated;
			goto out;
		}

		step++;
		bytes = sb->s_blocksize * allocated - offset;
		if (bytes > count)
			bytes = count;

		kmem = nova_get_block(inode->i_sb,
			nova_get_block_off(sb, blocknr,	pi->i_blk_type));

		if (offset || ((offset + bytes) & (PAGE_SIZE - 1)) != 0)
			nova_handle_head_tail_blocks(sb, pi, inode, pos, bytes,
								kmem);

		/* Now copy from user buf */
//		nova_dbg("Write: %p\n", kmem);
		NOVA_START_TIMING(memcpy_w_nvmm_t, memcpy_time);
		copied = bytes - memcpy_to_pmem_nocache(kmem + offset,
						buf, bytes);
		NOVA_END_TIMING(memcpy_w_nvmm_t, memcpy_time);

		entry_data.pgoff = cpu_to_le64(start_blk);
		entry_data.num_pages = cpu_to_le32(allocated);
		entry_data.invalid_pages = 0;
		entry_data.block = cpu_to_le64(nova_get_block_off(sb, blocknr,
							pi->i_blk_type));
		entry_data.mtime = cpu_to_le32(time);
		/* Set entry type after set block */
		nova_set_entry_type((void *)&entry_data, FILE_WRITE);

		if (pos + copied > inode->i_size)
			entry_data.size = cpu_to_le64(pos + copied);
		else
			entry_data.size = cpu_to_le64(inode->i_size);

		curr_entry = nova_append_file_write_entry(sb, pi, inode,
							&entry_data, temp_tail);
		if (curr_entry == 0) {
			nova_dbg("%s: append inode entry failed\n", __func__);
			ret = -ENOSPC;
			goto out;
		}

		nova_dbgv("Write: %p, %lu\n", kmem, copied);
		if (copied > 0) {
			status = copied;
			written += copied;
			pos += copied;
			buf += copied;
			count -= copied;
			num_blocks -= allocated;
		}
		if (unlikely(copied != bytes)) {
			nova_dbg("%s ERROR!: %p, bytes %lu, copied %lu\n",
				__func__, kmem, bytes, copied);
			if (status >= 0)
				status = -EFAULT;
		}
		if (status < 0)
			break;

		if (begin_tail == 0)
			begin_tail = curr_entry;
		temp_tail = curr_entry + sizeof(struct nova_file_write_entry);
	}

	nova_memunlock_inode(sb, pi);
	data_bits = blk_type_to_shift[pi->i_blk_type];
	le64_add_cpu(&pi->i_blocks,
			(total_blocks << (data_bits - sb->s_blocksize_bits)));
	nova_memlock_inode(sb, pi);

	nova_update_tail(pi, temp_tail);

	/* Free the overlap blocks after the write is committed */
	ret = nova_reassign_file_tree(sb, pi, sih, begin_tail);
	if (ret)
		goto out;

	inode->i_blocks = le64_to_cpu(pi->i_blocks);

	ret = written;
	NOVA_STATS_ADD(write_breaks, step);
	nova_dbgv("blocks: %lu, %llu\n", inode->i_blocks, pi->i_blocks);

	*ppos = pos;
	if (pos > inode->i_size) {
		i_size_write(inode, pos);
		sih->i_size = pos;
	}

out:
	if (ret < 0)
		nova_cleanup_incomplete_write(sb, pi, sih, blocknr, allocated,
						begin_tail, temp_tail);

	if (need_mutex)
		mutex_unlock(&inode->i_mutex);
	sb_end_write(inode->i_sb);
	NOVA_END_TIMING(cow_write_t, cow_write_time);
	NOVA_STATS_ADD(cow_write_bytes, written);
	return ret;
}

#ifdef JOURNAL_WRITE
int nova_reassign_file_tree_nofree(struct super_block *sb,
	struct nova_inode *pi, struct nova_inode_info_header *sih,
	u64 begin_tail, struct nova_file_write_entry *old_entry1,
	struct nova_file_write_entry *old_entry2)
{
	struct nova_file_write_entry *entry_data;
	u64 curr_p = begin_tail;
	size_t entry_size = sizeof(struct nova_file_write_entry);
	
	//first write entry	
	if (curr_p != pi->log_tail) {
		if (is_last_entry(curr_p, entry_size))
			curr_p = next_log_page(sb, curr_p);

		if (curr_p == 0) {
			nova_err(sb, "%s: File inode %llu log is NULL!\n",
				__func__, pi->nova_ino);
			return -EINVAL;
		}

		entry_data = (struct nova_file_write_entry *)
					nova_get_block(sb, curr_p);

		if (nova_get_entry_type(entry_data) != FILE_WRITE) {
			nova_dbg("%s: entry type is not write? %d\n",
				__func__, nova_get_entry_type(entry_data));
			curr_p += entry_size;
			//continue;
		}

		nova_update_write_entry(sb, pi, sih, entry_data, old_entry1);
		curr_p += entry_size;
	}

	//2nd write entry
	if (curr_p != pi->log_tail) {
		if (old_entry2 == NULL) {
			nova_dbg("%s: need to update 2 write entry but 2nd old entry is NULL! \n", __func__);
		}
			
		if (is_last_entry(curr_p, entry_size))
			curr_p = next_log_page(sb, curr_p);

		if (curr_p == 0) {
			nova_err(sb, "%s: File inode %llu log is NULL!\n",
				__func__, pi->nova_ino);
			return -EINVAL;
		}

		entry_data = (struct nova_file_write_entry *)
					nova_get_block(sb, curr_p);

		if (nova_get_entry_type(entry_data) != FILE_WRITE) {
			nova_dbg("%s: entry type is not write? %d\n",
				__func__, nova_get_entry_type(entry_data));
			curr_p += entry_size;
			//continue;
		}

		nova_update_write_entry(sb, pi, sih, entry_data, old_entry2);
		curr_p += entry_size;
	}
	
	if(curr_p != pi->log_tail) {
		nova_dbg("%s have more than 2 write entry \n",__func__);
	}

	return 0;
}

#if 0
int nova_reassign_file_tree_nofree(struct super_block *sb,
	struct nova_inode *pi, struct nova_inode_info_header *sih,
	u64 begin_tail)
{
	struct nova_file_write_entry *entry_data;
	u64 curr_p = begin_tail;
	size_t entry_size = sizeof(struct nova_file_write_entry);

	while (curr_p != pi->log_tail) {
		if (is_last_entry(curr_p, entry_size))
			curr_p = next_log_page(sb, curr_p);

		if (curr_p == 0) {
			nova_err(sb, "%s: File inode %llu log is NULL!\n",
				__func__, pi->nova_ino);
			return -EINVAL;
		}

		entry_data = (struct nova_file_write_entry *)
					nova_get_block(sb, curr_p);

		if (nova_get_entry_type(entry_data) != FILE_WRITE) {
			nova_dbg("%s: entry type is not write? %d\n",
				__func__, nova_get_entry_type(entry_data));
			curr_p += entry_size;
			continue;
		}

		nova_assign_write_entry(sb, pi, sih, entry_data, false);
		curr_p += entry_size;
	}

	return 0;
}
#endif
bool nova_dax_journal_write_check(struct file *filp, const char __user *buf,
	size_t len, loff_t *ppos)
{

	struct address_space *mapping = filp->f_mapping;
	struct inode    *inode = mapping->host;
	struct super_block *sb = inode->i_sb;
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;

	loff_t pgoff_start=(*ppos)>> sb->s_blocksize_bits;//block start
	loff_t pgoff_end=(*ppos+len)>> sb->s_blocksize_bits;
	loff_t pgoff_file=inode->i_size>>sb->s_blocksize_bits;
	void **entry1=NULL;
	//void **entry2=NULL;

	//pgoff_start=(*ppos)>> sb->s_blocksize_bits;//block start
	//pgoff_end=(*ppos+len)>> sb->s_blocksize_bits;

	if (WH_DEBUG == 2) {
		nova_dbg("%s is called\n",__func__);
	}
		
	if (inode->i_size == 0 || (pgoff_end > pgoff_file)) {
		if (WH_DEBUG == 2) {
			nova_dbg("%s extend file size; return false\n",__func__);
			nova_dbg("%s :pgoff_end %llu ; pgoff_file %llu \n",__func__, pgoff_end, pgoff_file);
		}
		return false;
	}
	
#ifdef DYNAMIC_JOURNAL
	write_len_sum += len;
	if(write_count == 1024) {
		write_len_average = write_len_sum>>10 > NOVA_DATA_JOURNAL_SIZE? NOVA_DATA_JOURNAL_SIZE: write_len_sum>>10;
		write_len_sum = 0;
		write_count = 0;
	} else {
		write_count++;
	}

	if (len>write_len_average) {
#else
	if (len>NOVA_DATA_JOURNAL_SIZE) {
#endif
		if (WH_DEBUG == 2) {
			nova_dbg("%s len is too large; return false\n",__func__);
		}
		return false;
	}
//	if (filp->f_flags & O_APPEND){
//		if (WH_DEBUG == 2) {
//			nova_dbg("%s append write; return false\n",__func__);
//		}
//		return false;
//	}

	//check whether there is a write_entry for page to be written	
	entry1 = radix_tree_lookup_slot(&sih->tree, pgoff_start);
	if (entry1 == NULL) {
		if (WH_DEBUG == 2) {
			nova_dbg("%s entry not found; return false\n",__func__);
		}
		return false;
	}
	
	if (pgoff_start != pgoff_end) {
		return false;
	//	entry2 = radix_tree_lookup_slot(&sih->tree, pgoff_end);	
	//	if(entry2 == NULL){
	//		if (WH_DEBUG == 2) {
	//			nova_dbg("%s entry 2 not found; return false\n",__func__);
	//		}
	//		return false;
	//	}
	}
	//if there is one hole, give up journal write.
	
	//both or current page is available for inplace update
	if (WH_DEBUG == 2) {
		nova_dbg("%s return true\n",__func__);
	}
	return true;

}

ssize_t nova_dax_journal_write(struct file *filp,
	const char __user *buf,	size_t len, loff_t *ppos, bool need_mutex)
{
	struct address_space *mapping = filp->f_mapping;
	struct inode    *inode = mapping->host;
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	struct super_block *sb = inode->i_sb;
	//struct nova_sb_info *sbi = NOVA_SB(sb);
	struct nova_inode *pi;
	void **pentry1, **pentry2;
	struct nova_file_write_entry *old_entry1, *old_entry2; //*update_entry1, *update_entry2;
	struct nova_file_write_entry entry_data1, entry_data2;//mtime changes, need to replace write entry

	struct nova_data_journal_entry journal_entry_data1,journal_entry_data2;
	struct ptr_pair *journal_ptr;
	
	loff_t pos;
	size_t count1, count2, offset;
	size_t copied = 0;
	long status = 0;

	loff_t pgoff_start = (*ppos) >> sb->s_blocksize_bits;//block start
	loff_t pgoff_end = (*ppos + len) >> sb->s_blocksize_bits;
	//loff_t pgoff_max;

	timing_t journal_write_time;
	u64 temp_tail = 0, begin_tail = 0;
	u64 curr_entry;
	u64 temp;
	u32 time;
	int cpuid;
	int ret;
	unsigned long nvmm;
	void *ptr1, *ptr2;	
	
	if (WH_DEBUG == 1) {
		nova_dbg("%s: pgoff_start=%llu, pgoff_end=%llu, len=%zu", __func__,
				pgoff_start, pgoff_end, len);
	}
	
	
	if (len == 0)
		return 0;

	//Since it is update in place, it can be mmaped
	//do not check mapping
	NOVA_START_TIMING(journal_write_t, journal_write_time);
	
	sb_start_write(inode->i_sb);
	if (need_mutex)
		mutex_lock(&inode->i_mutex);
	
	if (!access_ok(VERIFY_READ, buf, len)) {
		ret = -EFAULT;
		goto out;
	}	

	pos = *ppos;
	
	if (filp->f_flags & O_APPEND)
		pos=i_size_read(inode);
	//O_APPEND is checked in nova_dax_journal_write_check()

	
	pi = nova_get_inode(sb, inode);
	//offset = pos & (sb->s_blocksize - 1);//relative offset in block
	offset = pos & (nova_inode_blk_size(pi) - 1);
	
	if (WH_DEBUG == 1) {
		nova_dbg("%s: pos=%llu, offset=%zu", __func__, pos, offset);
	}
	
	ret = file_remove_privs(filp);
	if (ret) {
		goto out;
	}
	
	inode->i_ctime = inode->i_mtime = CURRENT_TIME_SEC;
	time = CURRENT_TIME_SEC.tv_sec;
	
	nova_dbgv("%s: inode %lu, offset %lld, count %lu\n",
			__func__, inode->i_ino,	pos, len);
	
	//start to journal
	//get first write entry
	pentry1 = radix_tree_lookup_slot(&sih->tree, pgoff_start);
		
	//pgoff_max=old_entry1->pgoff+old_entry1->num_pages;

	//May need another write entry
	if (pgoff_start != pgoff_end) {
		if (WH_DEBUG == 1) {
			nova_dbg("%s: \n pgoff_start %lld != pgoff_end %lld \n",
					__func__, pgoff_start, pgoff_end);
		}
		pentry2 = radix_tree_lookup_slot(&sih->tree, pgoff_end);
		if(pentry2 == pentry1){
			pentry2 = NULL;
			if (WH_DEBUG == 1) {
				nova_dbg("two block in the same write entry \n");
			}
		}	
	} else {
		pentry2 = NULL;
	}

	old_entry1 = radix_tree_deref_slot(pentry1);
	
	if (WH_DEBUG == 1) {
		nova_dbg("old entry 1 with pgoff=%llu num_pages=%u block=%llu",
				old_entry1->pgoff, old_entry1->num_pages, old_entry1->block);
	}

	if (pentry2) {
		old_entry2 = radix_tree_deref_slot(pentry2);
	} else {
		old_entry2 = NULL;
	}
	//if entry2==NULL, we only need one journal entry
	//build entry_data1
	//get block
	if (WH_DEBUG == 1) {
		nova_dbg("%s: get_nvmm called\n",__func__);
	}
	nvmm = get_nvmm(sb, sih, old_entry1, pgoff_start);
	
	/*start block byte*/
	ptr1 = nova_get_block(sb, (nvmm << PAGE_SHIFT)+offset);

	journal_entry_data1.addr = (u64)(nvmm << PAGE_SHIFT)+offset;
	count1 = journal_entry_data1.length = old_entry2 != NULL? nova_inode_blk_size(pi) - offset : len;
	//nova_dbg("write entry 1 count1=%zu \n",count1);
	
	//copy data; Or use __copy_from_user_inatomic_nocache?
	copy_to_user(journal_entry_data1.data, ptr1, count1);
	
	//fill rest with zero
	//memset(journal_entry_data1->data + count1, 0, NOVA_DATA_JOURNAL_SIZE - count1);
	
	//build and append journal entry2
	if (old_entry2) {
		nvmm = get_nvmm(sb, sih, old_entry2, pgoff_end);
		/*start block byte*/
		ptr2 = nova_get_block(sb, nvmm << PAGE_SHIFT);		
		
		journal_entry_data2.addr=(u64)nvmm << PAGE_SHIFT;
		count2 = journal_entry_data2.length=((pos + len) & (nova_inode_blk_size(pi) - 1));
		//nova_dbg("write entry 2 count1=%zu \n",count2);

		copy_to_user(journal_entry_data2.data, ptr2, count2);
		//memset(journal_entry_data2->data + count2, 0, NOVA_DATA_JOURNAL_SIZE - count2);
	}

	//Journal old data
	cpuid = smp_processor_id();
	journal_ptr = nova_get_data_journal_pointers(sb,cpuid);
	temp = journal_ptr->journal_tail;	

	nova_append_data_journal_entry(sb, &journal_entry_data1, temp);
	temp = next_data_journal(temp);
	
	if (old_entry2) {
		nova_append_data_journal_entry(sb, &journal_entry_data2, temp); 
		temp = next_data_journal(temp);
	}
	journal_ptr->journal_tail = temp;

	//Memcpy write
	while (copied != count1) {
		copied = count1 - memcpy_to_pmem_nocache(ptr1, buf + copied, count1 - copied);
		if (WH_DEBUG == 1) {
			nova_dbg("Memcpy count1 =%zu copied = %zu", count1, copied);
		}
	}

	//Build write entry 1
	entry_data1.pgoff = cpu_to_le64(old_entry1->pgoff);
	entry_data1.num_pages = cpu_to_le32(old_entry1->num_pages);
	entry_data1.invalid_pages = cpu_to_le32(old_entry1->invalid_pages);
	entry_data1.block = cpu_to_le64(old_entry1->block);
	//entry_data1.size = cpu_to_le64(old_entry1->size);
	
	if (pos + copied > inode->i_size)
		entry_data1.size = cpu_to_le64(pos + copied);
	else
		entry_data1.size = cpu_to_le64(inode->i_size);
	
	//Only change modifed time
	entry_data1.mtime = cpu_to_le32(time);
	nova_set_entry_type((void *)&entry_data1, FILE_WRITE);
	if (WH_DEBUG == 1) {
		nova_dbg("%s: \n build write entry 1: pgoff=%llu, num_pages=%u, block=%llu \n",
				__func__, entry_data1.pgoff, entry_data1.num_pages, entry_data1.block);
	}
	
	temp_tail = pi->log_tail;
	curr_entry = nova_append_file_write_entry(sb, pi, inode,
							&entry_data1, temp_tail);	

	if (!curr_entry) {
		nova_dbg("%s: append inode entry failed\n", __func__);
		ret = -ENOSPC;
		goto undo;
	}

	if (begin_tail == 0)
		begin_tail = curr_entry;
	
	temp_tail = curr_entry + sizeof(struct nova_file_write_entry);
	
	if (old_entry2) {
		while (copied != len){
			copied += count2 - memcpy_to_pmem_nocache(ptr2, buf+copied, len-copied);
			if (WH_DEBUG == 1) {
				nova_dbg("Memcpy entry 2 copied = %zu remain bytes = %zu", copied, len - copied);
			}
		}
		//Build write entry 2
		entry_data2.pgoff = cpu_to_le64(old_entry2->pgoff);
		entry_data2.num_pages = cpu_to_le32(old_entry2->num_pages);
		entry_data2.invalid_pages = cpu_to_le32(old_entry2->invalid_pages);
		entry_data2.block = cpu_to_le64(old_entry2->block);
		//Only change modifed time
		entry_data2.mtime = cpu_to_le32(time);
		entry_data2.size=cpu_to_le64(old_entry2->size);
	
		if (WH_DEBUG == 1) {
			nova_dbg("finish write entry 2 with pgoff=%llu num_pages=%u block=%llu",
					entry_data2.pgoff, entry_data2.num_pages, entry_data2.block);
		}
		curr_entry = nova_append_file_write_entry(sb, pi, inode,
							&entry_data2, temp_tail);			

		if (!curr_entry) {
			nova_dbg("%s: append inode entry failed\n", __func__);
			ret = -ENOSPC;
			goto undo;
		}
		temp_tail = curr_entry + sizeof(struct nova_file_write_entry);
	}

	if (copied > 0) {
		pos+=copied;
	}

	if (unlikely(copied != len)) {
		nova_dbg("%s ERROR!: %p, bytes %lu, copied %lu\n",
		__func__, ptr1 , len, copied);
		if (status >= 0)
			status = -EFAULT;
	}

	
	
	nova_update_tail(pi, temp_tail);
	
	//Update radix tree without free data block
	//update_entry1 = (struct nova_file_write_entry *)nova_get_block(sb, (u64)&entry_data1);
	
	ret = nova_reassign_file_tree_nofree(sb, pi, sih, begin_tail, old_entry1, old_entry2);
	
	//ret = nova_assign_write_entry(sb, pi, sih, update_entry1, false);
	//if (old_entry2){
	//	update_entry2 = (struct nova_file_write_entry *)
	//				nova_get_block(sb, (u64)&entry_data2);
		//ret += nova_assign_write_entry(sb, pi, sih, &entry_data2, false);
	//}
	
	if (ret) {
		nova_dbg("%s: ERROR %d\n", __func__, ret);
		goto undo;
	}

	ret=copied;	
	
	*ppos = pos;
	if (pos > inode->i_size) {
		i_size_write(inode, pos);
		sih->i_size = pos;
	}
	//Commit write change journal head to tail
	nova_commit_data_transaction(sb, temp, cpuid);
	
undo:
	//Append file entry failure undo 
	if (ret<0){
		if (old_entry2){
			nova_recover_data_journal(sb, journal_ptr, 2);	
		} else {
			nova_recover_data_journal(sb, journal_ptr, 1);
		}
	}

out:
	if (need_mutex)
		mutex_unlock(&inode->i_mutex);
	sb_end_write(inode->i_sb);
	NOVA_END_TIMING(journal_write_t, journal_write_time);	
	NOVA_STATS_ADD(journal_write_bytes, copied);
	if (WH_DEBUG == 1) {
		nova_dbg("%s: return %d \n", __func__, ret);
	}
	return ret;

}


ssize_t nova_dax_file_write(struct file *filp, const char __user *buf,
	size_t len, loff_t *ppos)
{	
	if (WH_DEBUG == 2) {
		nova_dbg("%s called \n", __func__);
	}
	if (nova_dax_journal_write_check(filp, buf, len, ppos)){// start journal write
		return nova_dax_journal_write(filp, buf, len, ppos, true);
	} else {
		return nova_cow_file_write(filp, buf, len, ppos, true);
	}
}

#endif /*JOURNAL_WRITE*/ 

#ifndef JOURNAL_WRITE

ssize_t nova_dax_file_write(struct file *filp, const char __user *buf,
	size_t len, loff_t *ppos)
{	
	return nova_cow_file_write(filp, buf, len, ppos, true);
	
}


#endif

#if 0
ssize_t nova_dax_file_write(struct file *filp, const char __user *buf,
	size_t len, loff_t *ppos)
{
	return nova_cow_file_write(filp, buf, len, ppos, true);
}
#endif

/*
 * return > 0, # of blocks mapped or allocated.
 * return = 0, if plain lookup failed.
 * return < 0, error case.
 */
static int nova_dax_get_blocks(struct inode *inode, sector_t iblock,
	unsigned long max_blocks, struct buffer_head *bh, int create)
{
	struct super_block *sb = inode->i_sb;
	struct nova_inode *pi;
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	struct nova_file_write_entry *entry = NULL;
	struct nova_file_write_entry entry_data;
	u64 temp_tail = 0;
	u64 curr_entry;
	u32 time;
	unsigned int data_bits;
	unsigned long nvmm = 0;
	unsigned long next_pgoff;
	unsigned long blocknr = 0;
	int num_blocks = 0;
	int allocated = 0;
	int ret = 0;

	if (max_blocks == 0)
		return 0;

	nova_dbgv("%s: pgoff %lu, num %lu, create %d\n",
				__func__, iblock, max_blocks, create);

	entry = nova_get_write_entry(sb, si, iblock);
	if (entry) {
		/* Find contiguous blocks */
		if (entry->invalid_pages == 0)
			num_blocks = entry->num_pages - (iblock - entry->pgoff);
		else
			num_blocks = 1;

		if (num_blocks > max_blocks)
			num_blocks = max_blocks;

		nvmm = get_nvmm(sb, sih, entry, iblock);
		clear_buffer_new(bh);
		nova_dbgv("%s: pgoff %lu, block %lu\n", __func__, iblock, nvmm);
		goto out;
	}

	if (create == 0)
		return 0;

	pi = nova_get_inode(sb, inode);
	num_blocks = max_blocks;
	inode->i_ctime = inode->i_mtime = CURRENT_TIME_SEC;
	time = CURRENT_TIME_SEC.tv_sec;

	/* Fill the hole */
	entry = nova_find_next_entry(sb, sih, iblock);
	if (entry) {
		next_pgoff = entry->pgoff;
		if (next_pgoff <= iblock) {
			BUG();
			ret = -EINVAL;
			goto out;
		}

		num_blocks = next_pgoff - iblock;
		if (num_blocks > max_blocks)
			num_blocks = max_blocks;
	}

	/* Return initialized blocks to the user */
	allocated = nova_new_data_blocks(sb, pi, &blocknr, num_blocks,
						iblock, 1, 1);
	if (allocated <= 0) {
		nova_dbg("%s alloc blocks failed %d\n", __func__,
							allocated);
		ret = allocated;
		goto out;
	}

	num_blocks = allocated;
	entry_data.pgoff = cpu_to_le64(iblock);
	entry_data.num_pages = cpu_to_le32(num_blocks);
	entry_data.invalid_pages = 0;
	entry_data.block = cpu_to_le64(nova_get_block_off(sb, blocknr,
							pi->i_blk_type));
	/* Set entry type after set block */
	nova_set_entry_type((void *)&entry_data, FILE_WRITE);
	entry_data.mtime = cpu_to_le32(time);

	/* Do not extend file size */
	entry_data.size = cpu_to_le64(inode->i_size);

	curr_entry = nova_append_file_write_entry(sb, pi, inode,
						&entry_data, pi->log_tail);
	if (curr_entry == 0) {
		nova_dbg("%s: append inode entry failed\n", __func__);
		ret = -ENOSPC;
		goto out;
	}

	nvmm = blocknr;
	data_bits = blk_type_to_shift[pi->i_blk_type];
	le64_add_cpu(&pi->i_blocks,
			(num_blocks << (data_bits - sb->s_blocksize_bits)));

	temp_tail = curr_entry + sizeof(struct nova_file_write_entry);
	nova_update_tail(pi, temp_tail);

	ret = nova_reassign_file_tree(sb, pi, sih, curr_entry);
	if (ret)
		goto out;

	inode->i_blocks = le64_to_cpu(pi->i_blocks);

//	set_buffer_new(bh);

out:
	if (ret < 0) {
		nova_cleanup_incomplete_write(sb, pi, sih, blocknr, allocated,
						0, temp_tail);
		return ret;
	}

	map_bh(bh, inode->i_sb, nvmm);
	if (num_blocks > 1)
		bh->b_size = sb->s_blocksize * num_blocks;

	return num_blocks;
}

int nova_dax_get_block(struct inode *inode, sector_t iblock,
	struct buffer_head *bh, int create)
{
	unsigned long max_blocks = bh->b_size >> inode->i_blkbits;
	int ret;
	timing_t gb_time;

	NOVA_START_TIMING(dax_get_block_t, gb_time);

	ret = nova_dax_get_blocks(inode, iblock, max_blocks, bh, create);
	if (ret > 0) {
		bh->b_size = ret << inode->i_blkbits;
		ret = 0;
	}
	NOVA_END_TIMING(dax_get_block_t, gb_time);
	return ret;
}

#if 0
static ssize_t nova_flush_mmap_to_nvmm(struct super_block *sb,
	struct inode *inode, struct nova_inode *pi, loff_t pos,
	size_t count, void *kmem)
{
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	unsigned long start_blk;
	unsigned long cache_addr;
	u64 nvmm_block;
	void *nvmm_addr;
	loff_t offset;
	size_t bytes, copied;
	ssize_t written = 0;
	int status = 0;
	ssize_t ret;

	while (count) {
		start_blk = pos >> sb->s_blocksize_bits;
		offset = pos & (sb->s_blocksize - 1);
		bytes = sb->s_blocksize - offset;
		if (bytes > count)
			bytes = count;

		cache_addr = nova_get_cache_addr(sb, si, start_blk);
		if (cache_addr == 0) {
			nova_dbg("%s: ino %lu %lu mmap page %lu not found!\n",
					__func__, inode->i_ino, sih->ino, start_blk);
			nova_dbg("mmap pages %lu\n", sih->mmap_pages);
			ret = -EINVAL;
			goto out;
		}

		nvmm_block = MMAP_ADDR(cache_addr);
		nvmm_addr = nova_get_block(sb, nvmm_block);
		copied = bytes - memcpy_to_pmem_nocache(kmem + offset,
				nvmm_addr + offset, bytes);

		if (copied > 0) {
			status = copied;
			written += copied;
			pos += copied;
			count -= copied;
			kmem += offset + copied;
		}
		if (unlikely(copied != bytes)) {
			nova_dbg("%s ERROR!: %p, bytes %lu, copied %lu\n",
				__func__, kmem, bytes, copied);
			if (status >= 0)
				status = -EFAULT;
		}
		if (status < 0) {
			ret = status;
			goto out;
		}
	}
	ret = written;
out:
	return ret;
}

ssize_t nova_copy_to_nvmm(struct super_block *sb, struct inode *inode,
	struct nova_inode *pi, loff_t pos, size_t count, u64 *begin,
	u64 *end)
{
	struct nova_file_write_entry entry_data;
	unsigned long start_blk, num_blocks;
	unsigned long blocknr = 0;
	unsigned long total_blocks;
	unsigned int data_bits;
	int allocated = 0;
	u64 curr_entry;
	ssize_t written = 0;
	int ret;
	void *kmem;
	size_t bytes, copied;
	loff_t offset;
	int status = 0;
	u64 temp_tail = 0, begin_tail = 0;
	u32 time;
	timing_t memcpy_time, copy_to_nvmm_time;

	NOVA_START_TIMING(copy_to_nvmm_t, copy_to_nvmm_time);
	sb_start_write(inode->i_sb);

	offset = pos & (sb->s_blocksize - 1);
	num_blocks = ((count + offset - 1) >> sb->s_blocksize_bits) + 1;
	total_blocks = num_blocks;
	inode->i_ctime = inode->i_mtime = CURRENT_TIME_SEC;
	time = CURRENT_TIME_SEC.tv_sec;

	nova_dbgv("%s: ino %lu, block %llu, offset %lu, count %lu\n",
		__func__, inode->i_ino, pos >> sb->s_blocksize_bits,
		(unsigned long)offset, count);

	temp_tail = *end;
	while (num_blocks > 0) {
		offset = pos & (nova_inode_blk_size(pi) - 1);
		start_blk = pos >> sb->s_blocksize_bits;
		allocated = nova_new_data_blocks(sb, pi, &blocknr, num_blocks,
						start_blk, 0, 0);
		if (allocated <= 0) {
			nova_dbg("%s alloc blocks failed %d\n", __func__,
								allocated);
			ret = allocated;
			goto out;
		}

		bytes = sb->s_blocksize * allocated - offset;
		if (bytes > count)
			bytes = count;

		kmem = nova_get_block(inode->i_sb,
			nova_get_block_off(sb, blocknr,	pi->i_blk_type));

		if (offset || ((offset + bytes) & (PAGE_SIZE - 1)))
			nova_handle_head_tail_blocks(sb, pi, inode, pos,
							bytes, kmem);

		NOVA_START_TIMING(memcpy_w_wb_t, memcpy_time);
		copied = nova_flush_mmap_to_nvmm(sb, inode, pi, pos, bytes,
							kmem);
		NOVA_END_TIMING(memcpy_w_wb_t, memcpy_time);

		entry_data.pgoff = cpu_to_le64(start_blk);
		entry_data.num_pages = cpu_to_le32(allocated);
		entry_data.invalid_pages = 0;
		entry_data.block = cpu_to_le64(nova_get_block_off(sb, blocknr,
							pi->i_blk_type));
		/* FIXME: should we use the page cache write time? */
		entry_data.mtime = cpu_to_le32(time);
		/* Set entry type after set block */
		nova_set_entry_type((void *)&entry_data, FILE_WRITE);

		entry_data.size = cpu_to_le64(inode->i_size);

		curr_entry = nova_append_file_write_entry(sb, pi, inode,
						&entry_data, temp_tail);
		if (curr_entry == 0) {
			nova_dbg("%s: append inode entry failed\n", __func__);
			ret = -ENOSPC;
			goto out;
		}

		nova_dbgv("Write: %p, %ld\n", kmem, copied);
		if (copied > 0) {
			status = copied;
			written += copied;
			pos += copied;
			count -= copied;
			num_blocks -= allocated;
		}
		if (unlikely(copied != bytes)) {
			nova_dbg("%s ERROR!: %p, bytes %lu, copied %lu\n",
				__func__, kmem, bytes, copied);
			if (status >= 0)
				status = -EFAULT;
		}
		if (status < 0) {
			ret = status;
			goto out;
		}

		if (begin_tail == 0)
			begin_tail = curr_entry;
		temp_tail = curr_entry + sizeof(struct nova_file_write_entry);
	}

	nova_memunlock_inode(sb, pi);
	data_bits = blk_type_to_shift[pi->i_blk_type];
	le64_add_cpu(&pi->i_blocks,
			(total_blocks << (data_bits - sb->s_blocksize_bits)));
	nova_memlock_inode(sb, pi);
	inode->i_blocks = le64_to_cpu(pi->i_blocks);

	*begin = begin_tail;
	*end = temp_tail;

	ret = written;
out:
	if (ret < 0)
		nova_cleanup_incomplete_write(sb, pi, sih, blocknr, allocated,
						begin_tail, temp_tail);

	sb_end_write(inode->i_sb);
	NOVA_END_TIMING(copy_to_nvmm_t, copy_to_nvmm_time);
	return ret;
}

static int nova_get_nvmm_pfn(struct super_block *sb, struct nova_inode *pi,
	struct nova_inode_info *si, u64 nvmm, pgoff_t pgoff,
	vm_flags_t vm_flags, void **kmem, unsigned long *pfn)
{
	struct nova_inode_info_header *sih = &si->header;
	u64 mmap_block;
	unsigned long cache_addr = 0;
	unsigned long blocknr = 0;
	void *mmap_addr;
	void *nvmm_addr;
	int ret;

	cache_addr = nova_get_cache_addr(sb, si, pgoff);

	if (cache_addr) {
		mmap_block = MMAP_ADDR(cache_addr);
		mmap_addr = nova_get_block(sb, mmap_block);
	} else {
		ret = nova_new_data_blocks(sb, pi, &blocknr, 1,
						pgoff, 0, 1);

		if (ret <= 0) {
			nova_dbg("%s alloc blocks failed %d\n",
					__func__, ret);
			return ret;
		}

		mmap_block = blocknr << PAGE_SHIFT;
		mmap_addr = nova_get_block(sb, mmap_block);

		if (vm_flags & VM_WRITE)
			mmap_block |= MMAP_WRITE_BIT;

		nova_dbgv("%s: inode %lu, pgoff %lu, mmap block 0x%llx\n",
			__func__, sih->ino, pgoff, mmap_block);

		ret = radix_tree_insert(&sih->cache_tree, pgoff,
					(void *)mmap_block);
		if (ret) {
			nova_dbg("%s: ERROR %d\n", __func__, ret);
			return ret;
		}

		sih->mmap_pages++;
		if (nvmm) {
			/* Copy from NVMM to dram */
			nvmm_addr = nova_get_block(sb, nvmm);
			memcpy(mmap_addr, nvmm_addr, PAGE_SIZE);
		} else {
			memset(mmap_addr, 0, PAGE_SIZE);
		}
	}

	*kmem = mmap_addr;
	*pfn = nova_get_pfn(sb, mmap_block);

	return 0;
}

static int nova_get_mmap_addr(struct inode *inode, struct vm_area_struct *vma,
	pgoff_t pgoff, int create, void **kmem, unsigned long *pfn)
{
	struct super_block *sb = inode->i_sb;
	struct nova_inode_info *si = NOVA_I(inode);
	struct nova_inode_info_header *sih = &si->header;
	struct nova_inode *pi;
	u64 nvmm;
	vm_flags_t vm_flags = vma->vm_flags;
	int ret;

	pi = nova_get_inode(sb, inode);

	nvmm = nova_find_nvmm_block(sb, si, NULL, pgoff);

	ret = nova_get_nvmm_pfn(sb, pi, si, nvmm, pgoff, vm_flags,
						kmem, pfn);

	if (vm_flags & VM_WRITE) {
		if (pgoff < sih->low_dirty)
			sih->low_dirty = pgoff;
		if (pgoff > sih->high_dirty)
			sih->high_dirty = pgoff;
	}

	return ret;
}

/* OOM err return with dax file fault handlers doesn't mean anything.
 * It would just cause the OS to go an unnecessary killing spree !
 */
static int __nova_dax_file_fault(struct vm_area_struct *vma,
				  struct vm_fault *vmf)
{
	struct address_space *mapping = vma->vm_file->f_mapping;
	struct inode *inode = mapping->host;
	pgoff_t size;
	void *dax_mem;
	unsigned long dax_pfn = 0;
	int err;
	int ret = VM_FAULT_SIGBUS;

	mutex_lock(&inode->i_mutex);
	size = (i_size_read(inode) + PAGE_SIZE - 1) >> PAGE_SHIFT;
	if (vmf->pgoff >= size) {
		nova_dbg("[%s:%d] pgoff >= size(SIGBUS). vm_start(0x%lx),"
			" vm_end(0x%lx), pgoff(0x%lx), VA(%lx), size 0x%lx\n",
			__func__, __LINE__, vma->vm_start, vma->vm_end,
			vmf->pgoff, (unsigned long)vmf->virtual_address, size);
		goto out;
	}

	err = nova_get_mmap_addr(inode, vma, vmf->pgoff, 1,
						&dax_mem, &dax_pfn);
	if (unlikely(err)) {
		nova_dbg("[%s:%d] get_mmap_addr failed. vm_start(0x%lx),"
			" vm_end(0x%lx), pgoff(0x%lx), VA(%lx)\n",
			__func__, __LINE__, vma->vm_start, vma->vm_end,
			vmf->pgoff, (unsigned long)vmf->virtual_address);
		goto out;
	}

	nova_dbgv("%s flags: vma 0x%lx, vmf 0x%x\n",
			__func__, vma->vm_flags, vmf->flags);

	nova_dbgv("DAX mmap: inode %lu, vm_start(0x%lx), vm_end(0x%lx), "
			"pgoff(0x%lx), vma pgoff(0x%lx), "
			"VA(0x%lx)->PA(0x%lx)\n",
			inode->i_ino, vma->vm_start, vma->vm_end, vmf->pgoff,
			vma->vm_pgoff, (unsigned long)vmf->virtual_address,
			(unsigned long)dax_pfn << PAGE_SHIFT);

	if (dax_pfn == 0)
		goto out;

#if LINUX_VERSION_CODE >= KERNEL_VERSION(4, 5, 0)
	err = vm_insert_mixed(vma, (unsigned long)vmf->virtual_address,
		__pfn_to_pfn_t(dax_pfn, PFN_DEV));
#else
	err = vm_insert_mixed(vma, (unsigned long)vmf->virtual_address, dax_pfn);
#endif

	if (err == -ENOMEM)
		goto out;
	/*
	 * err == -EBUSY is fine, we've raced against another thread
	 * that faulted-in the same page
	 */
	if (err != -EBUSY)
		BUG_ON(err);

	ret = VM_FAULT_NOPAGE;

out:
	mutex_unlock(&inode->i_mutex);
	return ret;
}

static int nova_dax_file_fault(struct vm_area_struct *vma, struct vm_fault *vmf)
{
	int ret = 0;
	timing_t fault_time;

	NOVA_START_TIMING(mmap_fault_t, fault_time);
	ret = __nova_dax_file_fault(vma, vmf);
	NOVA_END_TIMING(mmap_fault_t, fault_time);
	return ret;
}
#endif

static int nova_dax_fault(struct vm_area_struct *vma, struct vm_fault *vmf)
{
	struct inode *inode = file_inode(vma->vm_file);
	int ret = 0;
	timing_t fault_time;

	NOVA_START_TIMING(mmap_fault_t, fault_time);

	mutex_lock(&inode->i_mutex);
	ret = dax_fault(vma, vmf, nova_dax_get_block, NULL);
	mutex_unlock(&inode->i_mutex);

	NOVA_END_TIMING(mmap_fault_t, fault_time);
	return ret;
}

static int nova_dax_pmd_fault(struct vm_area_struct *vma, unsigned long addr,
	pmd_t *pmd, unsigned int flags)
{
	struct inode *inode = file_inode(vma->vm_file);
	int ret = 0;
	timing_t fault_time;

	NOVA_START_TIMING(mmap_fault_t, fault_time);

	mutex_lock(&inode->i_mutex);
	ret = dax_pmd_fault(vma, addr, pmd, flags, nova_dax_get_block, NULL);
	mutex_unlock(&inode->i_mutex);

	NOVA_END_TIMING(mmap_fault_t, fault_time);
	return ret;
}

static int nova_dax_pfn_mkwrite(struct vm_area_struct *vma,
	struct vm_fault *vmf)
{
	struct inode *inode = file_inode(vma->vm_file);
	loff_t size;
	int ret = 0;
	timing_t fault_time;

	NOVA_START_TIMING(mmap_fault_t, fault_time);

	mutex_lock(&inode->i_mutex);
	size = (i_size_read(inode) + PAGE_SIZE - 1) >> PAGE_SHIFT;
	if (vmf->pgoff >= size)
		ret = VM_FAULT_SIGBUS;
	else
		ret = dax_pfn_mkwrite(vma, vmf);
	mutex_unlock(&inode->i_mutex);

	NOVA_END_TIMING(mmap_fault_t, fault_time);
	return ret;
}

static const struct vm_operations_struct nova_dax_vm_ops = {
	.fault	= nova_dax_fault,
	.pmd_fault = nova_dax_pmd_fault,
	.page_mkwrite = nova_dax_fault,
	.pfn_mkwrite = nova_dax_pfn_mkwrite,
};

int nova_dax_file_mmap(struct file *file, struct vm_area_struct *vma)
{
	file_accessed(file);

	vma->vm_flags |= VM_MIXEDMAP | VM_HUGEPAGE;

	vma->vm_ops = &nova_dax_vm_ops;
	nova_dbg_mmap4k("[%s:%d] MMAP 4KPAGE vm_start(0x%lx),"
			" vm_end(0x%lx), vm_flags(0x%lx), "
			"vm_page_prot(0x%lx)\n", __func__,
			__LINE__, vma->vm_start, vma->vm_end,
			vma->vm_flags, pgprot_val(vma->vm_page_prot));

	return 0;
}
