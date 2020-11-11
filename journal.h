/*
 * NOVA journal header
 *
 * Copyright 2015-2016 Regents of the University of California,
 * UCSD Non-Volatile Systems Lab, Andiry Xu <jix024@cs.ucsd.edu>
 * Copyright 2012-2013 Intel Corporation
 * Copyright 2009-2011 Marco Stornelli <marco.stornelli@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms and conditions of the GNU General Public License,
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 
 * 51 Franklin St - Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef __NOVA_JOURNAL_H__
#define __NOVA_JOURNAL_H__
#include <linux/slab.h>

#define NOVA_DATA_JOURNAL_ENTRY_SIZE 1024
#define NOVA_DATA_JOURNAL_SIZE (NOVA_DATA_JOURNAL_ENTRY_SIZE - 2*sizeof(u64))//Wang 256- 8B addr -8B length
//struct ptr_pair;
/* Lite journal */
struct nova_lite_journal_entry {
	/* The highest byte of addr is type */
	u64 addrs[4];
	u64 values[4];
};

struct nova_data_journal_entry {
	u64 addr;
	u64 length;//0-255
	u8 data[NOVA_DATA_JOURNAL_SIZE];
};

int nova_lite_journal_soft_init(struct super_block *sb);
int nova_lite_journal_hard_init(struct super_block *sb);
u64 nova_create_lite_transaction(struct super_block *sb,
	struct nova_lite_journal_entry *dram_entry1,
	struct nova_lite_journal_entry *dram_entry2,
	int entries, int cpu);
void nova_commit_lite_transaction(struct super_block *sb, u64 tail, int cpu);

//Wang
u64 next_data_journal(u64 curr_p);

void nova_recover_data_journal_entry(struct super_block *sb,
	u64 addr, u8 *data, size_t len);

void nova_undo_data_journal_entry(struct super_block *sb,
	struct nova_data_journal_entry *entry);

u64 nova_append_data_journal_entry(struct super_block *sb, 
	struct nova_data_journal_entry *data_entry, u64 tail);

void nova_commit_data_transaction(struct super_block *sb, u64 tail, int cpu);


int nova_data_journal_soft_init(struct super_block *sb);

int nova_data_journal_hard_init(struct super_block *sb);
#endif    /* __NOVA_JOURNAL_H__ */
