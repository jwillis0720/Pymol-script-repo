import sys
import pandas
sys.path.append('/Users/jordanwillis/Mounts/vax/schief_scripts/pymol')
import data_2b

# import sqlite3 as sql
# import util
# from pymol import cmd,stored 
# import sys, re
# from pymol import stored



####TESTING PURPOSES






# conn = sql.connect("/Users/jordanwillis/Mounts/vax/pg9_fit/step3_make_germline_reversions_with_v2loop/output/analysis_all.db3")
# c = conn.cursor()

# qt =  c.execute('''
# SELECT
# 	score_table.struct_id,
# 	score_table.tag,
# 	score_table.score_value,
# 	r.chain_id,
# 	r.residue_number,
# 	r.insertion_code,
# 	r.pdb_residue_number,
# 	res.name3,
# 	res.res_type,
# 	AVG(interface.dG) as "dG",
# 	AVG(interface.dSASA) as "dSASA",
# 	AVG(
# 		interface.relative_dSASA_fraction
# 	) as "relativeSASA",
# 	AVG(interface.dhSASA) as "dhSASA",
# 	count(interface.dG) AS "counts_at_interface"
# FROM
# 	(
# 		SELECT
# 			s.struct_id,
# 			s.tag,
# 			structure_scores.score_value
# 		FROM
# 			structures s
# 		INNER JOIN structure_scores ON structure_scores.struct_id = s.struct_id
# 		WHERE
# 			tag LIKE "%germ%"
# 		GROUP BY
# 			tag
# 		ORDER BY
# 			score_value ASC
# 		LIMIT 90
# 	) score_table
# INNER JOIN residue_pdb_identification r ON r.struct_id = score_table.struct_id
# INNER JOIN residues res ON res.resNum = r.residue_number
# AND res.struct_id = r.struct_id
# LEFT OUTER JOIN interface_residues interface ON interface.struct_id = r.struct_id
# AND interface.resNum = res.resNum
# GROUP BY
# 	r.residue_number''')

#query = [j for j in [i for i in qt]]
#descriptions = [i[0] for i in c.description]

#results = []
#for i in query:
# 	dict_i = {}
# 	for des,value in zip(descriptions,i):
# 		dict_i[des] = value
# 	results.append(dict_i)

#with open('somefile.txt','w') as f:
#	for i in results:
#		value = i['counts_at_interface']
#		print type(value)
#		if value is None:
##		f.write("{} {} {} {}\n".format(i['chain_id'],i['pdb_residue_number'],i['name3'],value))
		#sf.write("{} {} {} {}".format(i['chain_id'],i['pdb_residue_number'],i['name3'],i['counts_at_interface'))



# if __name__ == '__main__':
#     pdb_file = sys.argv[1]
#     b_dict = residue_data_extract(pdb_file)
#     for chain in sorted(b_dict):
#         for resi in sorted(b_dict[chain]):
#             b, resn = b_dict[chain][resi]
#             print "b-factors %s %s %s %s  new B='%s'" % (pdb_file, chain, resn, resi, b)
#     sys.exit()


#plug = []
#for re in results:
#	print plug.append(re['counts_at_interface'])


#for re in results:
#	cmd.alter("* and resi {}".format(re["resNum"]),"b={}".format(re['c']))
	#break