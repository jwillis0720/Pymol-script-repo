from multiprocessing import Process
import pymongo
import sys
import os
import time

#instance_string = "mongodb://schieflab:bryanisawesome@52.5.193.108"
#instance_string = "mongodb://root:Jiffylite42@vax.scripps.edu:9998"
instance_string = "mongodb://localhost"
instance_one = pymongo.MongoClient(instance_string)
_index = [('cdr3_len',pymongo.ASCENDING),
		  ('v_gene.fam',pymongo.ASCENDING),
		  ('v_gene.gene',pymongo.ASCENDING),
		  ('d_gene.fam',pymongo.ASCENDING),
		  ('d_gene.gene',pymongo.ASCENDING),
		  ('j_gene.gene',pymongo.ASCENDING)]
#_index = [('raw_input',pymongo.TEXT)]

def make_index(db,coll):
	print "sleeping _id {}".format(os.getpid())
	time.sleep(5)
	#print _index
	#spawned_instance = pymongo.MongoClient(instance_string)[db][coll]
	#spawned_instance.create_index(_index)

database_names = instance_one.database_names()
for i in database_names:
	print i
	#if i[0] != 'd':
	#	database_names.remove(i)


#database_names = []
#database_names.append('d5684_H275MBCXX')
#database_names.append('d6471_H275NBCXX')
dbs_col = {}
for database in database_names:
	dbs_col[database] = instance_one[database].collection_names()


# print dbs_col
# sys.exit()
# for database in dbs_col:
# 	for col in dbs_col[database]:
# 		spawned_instance = pymongo.MongoClient(instance_string)[database][col]
# 		print database, col
# 		print spawned_instance.index_information()
# 	#for col in database_names[database]:
# 	#	print col
# 		#pymongo.MongoClient(instance_string)


jobs = []
for db in dbs_col:
	for col in dbs_col[db]:
		p = Process(target=make_index,args=(db,col))
		#p.start()
		jobs.append(p)
job_counter = 0
jobs_running = []
while jobs:
	for i in jobs:
		if job_counter < 11:
			print "starting {}".format(i.name)
			i.start()
			job_counter += 1
			jobs.remove(i)
			jobs_running.append(i)
	for jobs_ in jobs_running:
		#jobs_.join()
		if jobs_.is_alive():
			continue
		else:
			print "print job {} is done, removing".format(jobs_.name)
			jobs_running.remove(jobs_)
			job_counter -= 1
    print "sleeping distributor"
	time.sleep(10)

#for job in jobs:
#		job.join()