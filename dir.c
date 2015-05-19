/*
 * BRIEF DESCRIPTION
 *
 * File operations for directories.
 *
 * Copyright 2012-2013 Intel Corporation
 * Copyright 2009-2011 Marco Stornelli <marco.stornelli@gmail.com>
 * Copyright 2003 Sony Corporation
 * Copyright 2003 Matsushita Electric Industrial Co., Ltd.
 * 2003-2004 (c) MontaVista Software, Inc. , Steve Longerbeam
 * This file is licensed under the terms of the GNU General Public
 * License version 2. This program is licensed "as is" without any
 * warranty of any kind, whether express or implied.
 */

#include <linux/fs.h>
#include <linux/pagemap.h>
#include "pmfs.h"

/*
 *	Parent is locked. We do not take lock for RB Tree operations.
 */

#define DT2IF(dt) (((dt) << 12) & S_IFMT)
#define IF2DT(sif) (((sif) & S_IFMT) >> 12)

/* ========================= RB Tree operations ============================= */

static int pmfs_rbtree_compare_find_by_name(struct super_block *sb,
	struct pmfs_dir_node *curr, const char *name, int namelen)
{
	struct pmfs_log_direntry *entry;
	int min_len;

	if (!curr || curr->nvmm == 0)
		BUG();

	entry = (struct pmfs_log_direntry *)pmfs_get_block(sb, curr->nvmm);
	min_len = namelen < entry->name_len ? namelen : entry->name_len;

	pmfs_dbg_verbose("%s: %s %s, entry @0x%lx\n", __func__,
				name, entry->name, curr->nvmm);
	if (strncmp(name, entry->name, min_len) < 0)
		return -1;
	if (strncmp(name, entry->name, min_len) > 0)
		return 1;

	if (namelen < entry->name_len)
		return -1;
	if (namelen > entry->name_len)
		return 1;
	return 0;
}

struct pmfs_dir_node *pmfs_find_dir_node_by_name(struct super_block *sb,
	struct pmfs_inode *pi, struct inode *inode, const char *name,
	unsigned long name_len)
{
	struct pmfs_inode_info *si = PMFS_I(inode);
	struct pmfs_inode_info_header *sih = &si->header;
	struct pmfs_dir_node *curr;
	struct rb_node *temp;
	int compVal;

	temp = sih->dir_tree.rb_node;
	while (temp) {
		curr = container_of(temp, struct pmfs_dir_node, node);
		compVal = pmfs_rbtree_compare_find_by_name(sb, curr,
							name, name_len);

		if (compVal == -1) {
			temp = temp->rb_left;
		} else if (compVal == 1) {
			temp = temp->rb_right;
		} else {
			return curr;
		}
	}

	return NULL;
}

static inline struct pmfs_dir_node *pmfs_find_dir_node(struct super_block *sb,
	struct pmfs_inode *pi, struct inode *inode, struct dentry *dentry)
{
	const char *name = dentry->d_name.name;
	int namelen = dentry->d_name.len;

	return pmfs_find_dir_node_by_name(sb, pi, inode, name, namelen);
}

int pmfs_insert_dir_node_by_name(struct super_block *sb, struct pmfs_inode *pi,
	struct pmfs_inode_info_header *sih, const char *name, int namelen,
	u64 dir_entry)
{
	struct pmfs_dir_node *curr, *new;
	struct rb_node **temp, *parent;
	int compVal;

	pmfs_dbg_verbose("%s: insert %s @ 0x%llx\n", __func__, name, dir_entry);

	temp = &(sih->dir_tree.rb_node);
	parent = NULL;

	while (*temp) {
		curr = container_of(*temp, struct pmfs_dir_node, node);
		compVal = pmfs_rbtree_compare_find_by_name(sb, curr,
							name, namelen);
		parent = *temp;

		if (compVal == -1) {
			temp = &((*temp)->rb_left);
		} else if (compVal == 1) {
			temp = &((*temp)->rb_right);
		} else {
			pmfs_dbg("%s: entry %s already exists\n",
				__func__, name);
			return -EINVAL;
		}
	}

	new = pmfs_alloc_dirnode(sb);
	if (!new)
		return -ENOMEM;

	new->nvmm = dir_entry;
	rb_link_node(&new->node, parent, temp);
	rb_insert_color(&new->node, &sih->dir_tree);
//	pmfs_print_dir_tree(sb, inode);

	return 0;
}

static inline int pmfs_insert_dir_node(struct super_block *sb,
	struct pmfs_inode *pi, struct inode *inode, struct dentry *dentry,
	u64 dir_entry)
{
	struct pmfs_inode_info *si = PMFS_I(inode);
	struct pmfs_inode_info_header *sih = &si->header;
	const char *name = dentry->d_name.name;
	int namelen = dentry->d_name.len;

	return pmfs_insert_dir_node_by_name(sb, pi, sih, name,
						namelen, dir_entry);
}

void pmfs_remove_dir_node_by_name(struct super_block *sb, struct pmfs_inode *pi,
	struct pmfs_inode_info_header *sih, const char *name, int namelen)
{
	struct pmfs_dir_node *curr;
	struct rb_node *temp;
	int compVal;

	temp = sih->dir_tree.rb_node;
	while (temp) {
		curr = container_of(temp, struct pmfs_dir_node, node);
		compVal = pmfs_rbtree_compare_find_by_name(sb, curr, name,
								namelen);

		if (compVal == -1) {
			temp = temp->rb_left;
		} else if (compVal == 1) {
			temp = temp->rb_right;
		} else {
			rb_erase(&curr->node, &sih->dir_tree);
			pmfs_free_dirnode(sb, curr);
			break;
		}
	}

	return;
}

static inline void pmfs_remove_dir_node(struct super_block *sb,
	struct pmfs_inode *pi, struct inode *inode, struct dentry *dentry)
{
	struct pmfs_inode_info *si = PMFS_I(inode);
	struct pmfs_inode_info_header *sih = &si->header;
	const char *name = dentry->d_name.name;
	int namelen = dentry->d_name.len;

	return pmfs_remove_dir_node_by_name(sb, pi, sih, name, namelen);
}

void pmfs_print_dir_tree(struct super_block *sb, struct inode *inode)
{
	struct pmfs_inode_info *si = PMFS_I(inode);
	struct pmfs_inode_info_header *sih = &si->header;
	struct pmfs_dir_node *curr;
	struct pmfs_log_direntry *entry;
	struct rb_node *temp;

	pmfs_dbg("%s: dir ino %lu\n", __func__, inode->i_ino);
	temp = rb_first(&sih->dir_tree);
	while (temp) {
		curr = container_of(temp, struct pmfs_dir_node, node);

		if (!curr || curr->nvmm == 0)
			BUG();

		entry = (struct pmfs_log_direntry *)
				pmfs_get_block(sb, curr->nvmm);
		pmfs_dbg("%.*s\n", entry->name_len, entry->name);
		temp = rb_next(temp);
	}

	return;
}

void pmfs_delete_dir_tree(struct super_block *sb, struct inode *inode)
{
	struct pmfs_inode_info *si = PMFS_I(inode);
	struct pmfs_inode_info_header *sih = &si->header;
	struct pmfs_dir_node *curr;
	struct rb_node *temp;

	temp = rb_first(&sih->dir_tree);
	while (temp) {
		curr = container_of(temp, struct pmfs_dir_node, node);
		temp = rb_next(temp);
		rb_erase(&curr->node, &sih->dir_tree);
		pmfs_free_dirnode(sb, curr);
	}
	return;
}

/* ========================= Entry operations ============================= */

static int pmfs_add_dirent_to_buf(struct dentry *dentry, struct inode *inode,
	struct pmfs_inode *pidir, int inc_link)
{
	struct inode *dir = dentry->d_parent->d_inode;
	struct super_block *sb = dir->i_sb;
	const char *name = dentry->d_name.name;
	int namelen = dentry->d_name.len;
	unsigned short loglen;
	u64 curr_entry;
	u64 ino;

	pmfs_dbg_verbose("%s: %s %d\n", __func__, name, namelen);
	if (inode) {
		ino = inode->i_ino;
	} else {
		ino = 0;
	}
	/*
	 * XXX shouldn't update any times until successful
	 * completion of syscall, but too many callers depend
	 * on this.
	 */
	dir->i_mtime = dir->i_ctime = CURRENT_TIME_SEC;
	/*dir->i_version++; */

//	pmfs_memunlock_inode(dir->i_sb, pidir);
//	pidir->i_mtime = cpu_to_le32(dir->i_mtime.tv_sec);
//	pidir->i_ctime = cpu_to_le32(dir->i_ctime.tv_sec);
//	pmfs_memlock_inode(dir->i_sb, pidir);

	loglen = PMFS_DIR_LOG_REC_LEN(namelen);
	curr_entry = pmfs_append_dir_inode_entry(sb, pidir,
					dir, ino, dentry, loglen, 0, inc_link);
	pmfs_insert_dir_node(sb, pidir, dir, dentry, curr_entry);
	/* FIXME: Flush all data before update log_tail */
	pidir->log_tail = curr_entry + loglen;

	return 0;
}

/* adds a directory entry pointing to the inode. assumes the inode has
 * already been logged for consistency
 */
int pmfs_add_entry(pmfs_transaction_t *trans, struct dentry *dentry,
		struct inode *inode, int inc_link)
{
	struct inode *dir = dentry->d_parent->d_inode;
	struct super_block *sb = dir->i_sb;
	int retval = -EINVAL;
	struct pmfs_inode *pidir;
	timing_t add_entry_time;

	pmfs_dbg_verbose("%s: dir %lu new inode %lu\n", __func__, dir->i_ino,
				inode->i_ino);
	PMFS_START_TIMING(add_entry_t, add_entry_time);
	if (!dentry->d_name.len)
		return -EINVAL;

	pidir = pmfs_get_inode(sb, dir->i_ino);

	retval = pmfs_add_dirent_to_buf(dentry, inode, pidir, inc_link);
	PMFS_END_TIMING(add_entry_t, add_entry_time);
	return retval;
}

/* removes a directory entry pointing to the inode. assumes the inode has
 * already been logged for consistency
 */
int pmfs_remove_entry(pmfs_transaction_t *trans, struct dentry *dentry,
		struct inode *inode, int dec_link)
{
	struct super_block *sb = inode->i_sb;
	struct inode *dir = dentry->d_parent->d_inode;
	struct pmfs_inode *pidir;
	struct qstr *entry = &dentry->d_name;
	unsigned short loglen;
	u64 curr_entry;
	timing_t remove_entry_time;

	PMFS_START_TIMING(remove_entry_t, remove_entry_time);

	if (!dentry->d_name.len)
		return -EINVAL;

	pidir = pmfs_get_inode(sb, dir->i_ino);
	loglen = PMFS_DIR_LOG_REC_LEN(entry->len);
	curr_entry = pmfs_append_dir_inode_entry(sb, pidir, dir,
					0, dentry, loglen, 0, dec_link);
	/* FIXME: Flush all data before update log_tail */
	pidir->log_tail = curr_entry + loglen;
	pmfs_remove_dir_node(sb, pidir, dir, dentry);

	PMFS_END_TIMING(remove_entry_t, remove_entry_time);
	return 0;
}

inline int pmfs_replay_add_entry(struct super_block *sb, struct pmfs_inode *pi,
	struct pmfs_inode_info_header *sih, struct pmfs_log_direntry *entry,
	u64 curr_p)
{
	if (!entry->name_len)
		return -EINVAL;

	return pmfs_insert_dir_node_by_name(sb, pi, sih, entry->name,
					entry->name_len, curr_p);
}

inline int pmfs_replay_remove_entry(struct super_block *sb,
	struct pmfs_inode *pi, struct pmfs_inode_info_header *sih,
	struct pmfs_log_direntry *entry)
{
	pmfs_remove_dir_node_by_name(sb, pi, sih, entry->name,
					entry->name_len);
	return 0;
}

void pmfs_rebuild_dir_time_and_size(struct super_block *sb,
	struct pmfs_inode *pi, struct pmfs_log_direntry *entry)
{
	pi->i_ctime = cpu_to_le32(entry->ctime);
	pi->i_mtime = cpu_to_le32(entry->mtime);
	pi->i_size = cpu_to_le64(entry->size);
	pi->i_links_count = entry->links_count;
}

int pmfs_rebuild_dir_inode_tree(struct super_block *sb, struct pmfs_inode *pi,
	struct pmfs_inode_info_header *sih, unsigned long ino,
	struct scan_bitmap *bm)
{
	struct pmfs_log_direntry *entry;
	struct pmfs_inode_log_page *curr_page;
	u64 curr_p = pi->log_head;
	u64 next;
	int ret;

	pmfs_dbg_verbose("Rebuild dir %lu tree\n", ino);
	sih->dir_tree = RB_ROOT;

	if (curr_p == 0) {
		pmfs_err(sb, "Dir %lu log is NULL!\n", ino);
		BUG();
	}

	if (bm) {
		BUG_ON(curr_p & (PAGE_SIZE - 1));
		set_bit(curr_p >> PAGE_SHIFT, bm->bitmap_4k);
	}
	sih->log_pages = 1;
	while (curr_p != pi->log_tail) {
		if (curr_p == 0) {
			pmfs_err(sb, "Dir %lu log is NULL!\n", ino);
			BUG();
		}

		if (is_last_dir_entry(sb, curr_p)) {
			sih->log_pages++;
			curr_p = next_log_page(sb, curr_p);
			if (bm) {
				BUG_ON(curr_p & (PAGE_SIZE - 1));
				set_bit(curr_p >> PAGE_SHIFT, bm->bitmap_4k);
			}
		}

		pmfs_dbg_verbose("curr_p: 0x%llx\n", curr_p);
		entry = (struct pmfs_log_direntry *)pmfs_get_block(sb, curr_p);
		pmfs_dbg_verbose("entry @%p, ino %llu, name %*.s, namelen %u, "
			"rec len %u\n", entry, entry->ino, entry->name_len,
			entry->name, entry->name_len, entry->de_len);

		if (entry->ino > 0) {
			/* A valid entry to add */
			ret = pmfs_replay_add_entry(sb, pi, sih,
							entry, curr_p);
		} else {
			/* Delete the entry */
			ret = pmfs_replay_remove_entry(sb, pi, sih, entry);
		}
		pmfs_rebuild_dir_time_and_size(sb, pi, entry);
		curr_p += entry->de_len;
		if (ret) {
			pmfs_err(sb, "%s ERROR %d\n", __func__, ret);
			break;
		}
	}

	/* Keep traversing until log ends */
	curr_p &= PAGE_MASK;
	curr_page = (struct pmfs_inode_log_page *)pmfs_get_block(sb, curr_p);
	while ((next = curr_page->page_tail.next_page) != 0) {
		sih->log_pages++;
		curr_p = next;
		if (bm) {
			BUG_ON(curr_p & (PAGE_SIZE - 1));
			set_bit(curr_p >> PAGE_SHIFT, bm->bitmap_4k);
		}
		curr_page = (struct pmfs_inode_log_page *)
			pmfs_get_block(sb, curr_p);
	}

	return 0;
}

static int pmfs_readdir(struct file *file, struct dir_context *ctx)
{
	struct inode *inode = file_inode(file);
	struct super_block *sb = inode->i_sb;
	struct pmfs_inode *pi, *pidir;
	struct pmfs_inode_info *si = PMFS_I(inode);
	struct pmfs_inode_info_header *sih = &si->header;
	struct pmfs_dir_node *curr;
	struct pmfs_log_direntry *entry;
	struct rb_node *temp;
	ino_t ino;
	timing_t readdir_time;

	PMFS_START_TIMING(readdir_t, readdir_time);
	pidir = pmfs_get_inode(sb, inode->i_ino);
	pmfs_dbg_verbose("%s: ino %llu, root 0x%llx, size %llu, pos %llu\n",
				__func__, (u64)inode->i_ino, pidir->root,
				pidir->i_size, ctx->pos);

	if (ctx->pos == 0) {
		temp = rb_first(&sih->dir_tree);
	} else if (ctx->pos == READDIR_END) {
		goto out;
	} else if (ctx->pos) {
		entry = (struct pmfs_log_direntry *)
				pmfs_get_block(sb, ctx->pos);
		pmfs_dbg_verbose("ctx: ino %llu, name %*.s, "
				"name_len %u, de_len %u\n",
				(u64)entry->ino, entry->name_len, entry->name,
				entry->name_len, entry->de_len);
		curr = pmfs_find_dir_node_by_name(sb, NULL, inode,
					entry->name, entry->name_len);
		temp = &curr->node;
	}

	while (temp) {
		curr = container_of(temp, struct pmfs_dir_node, node);

		if (!curr || curr->nvmm == 0)
			BUG();

		entry = (struct pmfs_log_direntry *)
				pmfs_get_block(sb, curr->nvmm);
		if (entry->ino) {
			ino = le64_to_cpu(entry->ino);
			pi = pmfs_get_inode(sb, ino);
			pmfs_dbg_verbose("ctx: ino %llu, name %*.s, "
					"name_len %u, de_len %u\n",
					(u64)ino, entry->name_len, entry->name,
					entry->name_len, entry->de_len);
			if (!dir_emit(ctx, entry->name, entry->name_len,
					ino, IF2DT(le16_to_cpu(pi->i_mode)))) {
				pmfs_dbg_verbose("Here: pos %llu\n", ctx->pos);
				ctx->pos = curr->nvmm;
				return 0;
			}
		}
		temp = rb_next(temp);
	}

	/*
	 * We have reach the end. To let readdir be aware of that, we assign
	 * a bogus end offset to ctx.
	 */
	ctx->pos = READDIR_END;
out:
	PMFS_END_TIMING(readdir_t, readdir_time);
	return 0;
}

const struct file_operations pmfs_dir_operations = {
	.read		= generic_read_dir,
	.iterate	= pmfs_readdir,
	.fsync		= noop_fsync,
	.unlocked_ioctl = pmfs_ioctl,
#ifdef CONFIG_COMPAT
	.compat_ioctl	= pmfs_compat_ioctl,
#endif
};
