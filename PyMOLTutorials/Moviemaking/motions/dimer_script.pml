load 3f9f.pdb
as cartoon
util.cbc

create chainA, c. A
create chainB, c. B
select AnearB, br. c. A within 5 of c. B
select BnearA, br. c. B within 5 of c. A
deselect

delete 3f9f

orient
zoom *, 20

config_mouse three_button_motions

from pymol import movie

movie.add_blank(6)

set scene_buttons

