import sqlite3
import sys
from .mtree import *


class SQL_IO:
	dbName = ''
	nSnaps = 0	


	def __init__(self, dbName, nSnaps):	
		self.dbName = dbName	
		self.nSnaps = nSnaps
		self.mydb = sqlite3.connect(self.dbName)


	def halo_table(self):
		'Create halo table for SQL database'

		cursor = self.mydb.cursor()
		createStr = """CREATE TABLE IF NOT EXISTS halo (
				haloID INT64,
				simuCode VARCHAR(10), 
				allNumPart INT ARRAY[""" + str(self.nSnaps) + """], 
				allHaloIDs INT64 ARRAY[""" + str(self.nSnaps) + """],
				X FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				Y FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				Z FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				VX FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				VY FLOAT ARRAY[""" + str(self.nSnaps) + """], 
				VZ FLOAT ARRAY[""" + str(self.nSnaps) + """]
			)"""

		cursor.execute(createStr)
		self.mydb.commit()


	def insert_tree(self, ID, simuCode, haloParts, haloIDs):
		strCheck = self.select_tree(ID)
		
		if strCheck != None:
			print('ID: %s was already found in SQL database: %s' % (ID, self.dbName) )
		else:
			cursor = self.mydb.cursor()
			haloPartsStr = ''
			haloIDsStr = ''
		
			for iPart in range (0, len(haloParts)):
				if iPart < len(haloParts) - 1:
					haloPartsStr += str(haloParts[iPart]) + ', '
					haloIDsStr += str(haloIDs[iPart]) + ', '
				else:	
					haloPartsStr += str(haloParts[iPart])
					haloIDsStr += str(haloIDs[iPart])

			insertStr = """INSERT INTO halo (haloID, allNumPart, allHaloIDs, simuCode) 
					VALUES ('"""+ ID +"""', '""" + simuCode + """',  
					'""" + haloPartsStr + """', '""" + haloIDsStr + """' );"""

			#print(insertStr)
			cursor.execute(insertStr)
			self.mydb.commit()

	
	def select_tree(self, ID):
		cursor = self.mydb.cursor()
		selectStr = "SELECT * FROM halo WHERE haloID = %s " % str(ID)
		cursor.execute(selectStr)

		return cursor.fetchone()


	def get_full_mtree(self, ID):
		mergerTree = MTree(self.nSnaps, ID)
		
		tree_line = self.select_tree(ID)
		these_ids = tree_line[1].split()		
		these_nps = tree_line[3].split()		
		ids = []; nps = []; nPt = 0

		for nPt in range(0, self.nSnaps):
			
			try:
				this_id = these_ids[nPt].replace(",", "")
				this_np = these_nps[nPt].replace(",", "")
			except:
				#break
				'This is too short'
	
			if this_np == '' or this_np == ' ':
				this_np = 0	

			if this_id == '' or this_id == ' ':
				this_id = ids[nPt-1]
	
			ids.append(this_id)
			nps.append(this_np)

		mergerTree.fill_mass_id(nps, ids)

		return mergerTree

	def close(self):
		self.mydb.close()



