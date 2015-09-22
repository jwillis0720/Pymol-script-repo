create wt_revert_positions, WT and chain H and resi 100K+100L+100P+100R
create dj_revert_positions, *_model and chain H and resi 100K+100L+100P+100R
set grid_slot, 1, wt_revert_positions
set grid_slot, 2, dj_revert_positions
util.cbao wt_revert_positions
util.cbao dj_revert_positions
show spheres, *_revert_positions