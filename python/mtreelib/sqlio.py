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
		cursor = self.mydb.cursor()
		'Create halo table for SQL database'
		createStr = """
			CREATE TABLE halo (
			haloID INT64,
			allNumPart INT ARRAY[""" + str(self.nSnaps) + """], 
			allHaloIds INT64 ARRAY[""" + str(self.nSnaps) + """] 
			)"""

		cursor.execute(createStr)
		self.mydb.commit()


	def insert_tree(self, ID, nPart):
		cursor = self.mydb.cursor()

		nPartStr = ''
		
		for iPart in range (0, len(nPart)):
			if iPart < len(nPart) - 1:
				nPartStr += str(nPart[iPart]) + ', '
			else:	
				nPartStr += str(nPart[iPart])
	
		insertStr = """INSERT INTO halo (haloID, allNumPart) 
		VALUES (' """ + str(ID) + """ ', ' """ + nPartStr + """ ' );"""

		cursor.execute(insertStr)
		self.mydb.commit()

	
	def select_tree(self, ID):
		cursor = self.mydb.cursor()
		selectStr = "SELECT * FROM halo WHERE haloID = '%s' " % str(ID)
		cursor.execute(selectStr)
		print(cursor.fetchone())


	def close(self):
		self.mydb.close()



