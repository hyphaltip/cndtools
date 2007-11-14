# Copyright (c) 2006
# Colin Dewey (University of Wisconsin-Madison)
# cdewey@biostat.wisc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import _mysql
from _mysql import MySQLError
import operator

def isTupleOrList(x):
    """Return whether x is a tuple or a list."""
    return isinstance(x, (tuple, list))

class MySQLRecord(dict):
    """An object representation of a record from a database query.

       The field names of the record are attributes of the object.
       An object of this class can also be used as a dictionary.
    """
    def __init__(self, data):
        dict.__init__(self, data)
        self.__dict__ = data

class MySQLField:
    """A field from a database query."""
    def __init__(self, name, dataType, length):
        self.name = name
        self.type = dataType
        self.length = length
    def __repr__(self):
        return repr(self.__dict__)

class MySQLException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(self.value)

class MySQLResultIterator:
    def __init__(self, result, fmt):
        self.result = result
        self.fmt = fmt
    def __iter__(self):
        return self
    def next(self):
        rec = self.result.getNextRec(self.fmt)
        if rec is None:
            raise StopIteration
        else:
            return rec

class MySQLResult:
    FMT_TUPLE = 0
    FMT_DICT = 1
    FMT_LIST = 2
    FMT_OBJECT = 3

    def __init__(self,
                 resultObj,
                 query = None,
                 affectedRecs = None,
                 insertID = None):
        self.resultObj = resultObj
        self.query = query
        self.affectedRecs = affectedRecs
        self.insertID = insertID

    def __iter__(self):
        return MySQLResultIterator(self, fmt=MySQLResult.FMT_TUPLE)
    def getTupleIterator(self):
        return MySQLResultIterator(self, fmt=MySQLResult.FMT_TUPLE)
    def getDictIterator(self):
        return MySQLResultIterator(self, fmt=MySQLResult.FMT_DICT)
    def getListIterator(self):
        return MySQLResultIterator(self, fmt=MySQLResult.FMT_LIST)
    def getObjectIterator(self):
        return MySQLResultIterator(self, fmt=MySQLResult.FMT_OBJECT)

    def getQuery(self):
        if self.query is None:
            raise MySQLException("Error accessing query string:\n"
                                 + "Query string was not specified")
        return self.query

    def getNumAffectedRecs(self):
        if self.affectedRecs is None:
            raise MySQLException("Error accessing affected recs:\n"
                                 + "Query type does not affect rows\n"
                                 + "Query:\n"
                                 + self.query)
        return self.affectedRecs

    def getInsertID(self):
        if self.insertID is None:
            raise MySQLException("Error accessing insert id:\n"
                                 + "Query did not generate an"
                                 + "AUTO_INCREMENT value\n"
                                 + "Query:\n"
                                 + self.query)
        return self.insertID

    def getNextRec(self, fmt=FMT_TUPLE):
        rows = self.resultObj.fetch_row(how=(fmt % 2))
        if len(rows) > 0:
            if fmt == MySQLResult.FMT_LIST:
                return list(rows[0])
            elif fmt == MySQLResult.FMT_OBJECT:
                return MySQLRecord(rows[0])
            else:
                return rows[0]
        else:
            return None
    def getNextTuple(self):
        return self.getNextRec(fmt=MySQLResult.FMT_TUPLE)
    def getNextDict(self):
        return self.getNextRec(fmt=MySQLResult.FMT_DICT)
    def getNextList(self):
        return self.getNextRec(fmt=MySQLResult.FMT_LIST)
    def getNextObject(self):
        return self.getNextRec(fmt=MySQLResult.FMT_OBJECT)

    def getAllRecs(self, fmt=FMT_TUPLE):
        rows = self.resultObj.fetch_row(maxrows=0, how=(fmt % 2))
        if fmt == MySQLResult.FMT_LIST:
            return map(list, rows)
        elif fmt == MySQLResult.FMT_OBJECT:
            return map(MySQLRecord, rows)
        else:
            return rows
    def getAllTuples(self):
        return self.getAllRecs(fmt=MySQLResult.FMT_TUPLE)
    def getAllDicts(self):
        return self.getAllRecs(fmt=MySQLResult.FMT_DICT)
    def getAllLists(self):
        return self.getAllRecs(fmt=MySQLResult.FMT_LIST)
    def getAllObjects(self):
        return selt.getAllRecs(fmt=MySQLResult.FMT_OBJECT)

    def getNumRecs(self):
        return self.resultObj.num_rows()
    def getNumFields(self):
        return self.resultObj.num_fields()

    def seekToRec(self, rec_num):
        return self.resultObj.data_seek(rec_num)
    def getRecNum(self):
        return self.resultObj.row_tell()

    def _getFieldInfo(self, field_num, info_num):
        return self.resultObj.describe()[field_num][info_num]
    def _getAllFieldInfo(self, info_num):
        return map(operator.getitem,
                   self.resultObj.describe(),
                   [info_num] * self.getNumFields())
    def getField(self, fieldNum):
        return MySQLField(self.getFieldName(fieldNum),
                          self.getFieldType(fieldNum),
                          self.getFieldLen(fieldNum))
    def getFieldName(self, field_num):
        return self._getFieldInfo(field_num, 0)
    def getFieldType(self, field_num):
        return self._getFieldInfo(field_num, 1)
    def getFieldLen(self, field_num):
        return self._getFieldInfo(field_num, 3)

    def getAllFields(self):
        return [getField(i) for i in xrange(self.getNumFields())]
    def getAllFieldNames(self):
        return self._getAllFieldInfo(0)
    def getAllFieldTypes(self):
        return self._getAllFieldInfo(1)
    def getAllFieldLengths(self):
        return self._getAllFieldInfo(3)

class MySQLDB:
    from MySQLdb.constants import FIELD_TYPE

    conv = {FIELD_TYPE.TINY: int,
            FIELD_TYPE.SHORT: int,
            FIELD_TYPE.LONG: long,
            FIELD_TYPE.FLOAT: float,
            FIELD_TYPE.DOUBLE: float,
            FIELD_TYPE.DECIMAL: float,
            FIELD_TYPE.LONGLONG: long,
            FIELD_TYPE.INT24: int,
            FIELD_TYPE.YEAR: int
            }
    _current_db = None

    def __init__(self, host='', user='', passwd='', db=''):
        self.conn = _mysql.connect(host=host,
                                   user=user,
                                   passwd=passwd,
                                   db=db,
                                   conv=self.conv)
        self._current_db = db
    def stringLiteral(self, string):
        if string is None:
            return "NULL"
        else:
            return self.conn.string_literal(string)
    def setDB(self, db):
        self.conn.select_db(db)
        self._current_db = db
    def getDB(self):
        return self._current_db

    def makeTablename(self, db, table):
        return "%s.%s" % (db, table)

    def _query(self, query):
        try:
            self.conn.query(query)
        except MySQLError, e:
            raise MySQLException, "MySQL Error: " + str(e) + "\n" + \
                  "while attempting to execute query:\n" + query
        resultObj = self.conn.store_result()
        return MySQLResult(resultObj, query)

    def query(self, query):
        try:
            self.conn.query(query)
        except MySQLError, e:
            raise MySQLException, "MySQL Error: " + str(e) + "\n" + \
                  "while attempting to execute query:\n" + query        
        resultObj = self.conn.store_result()
        return MySQLResult(resultObj,
                           query,
                           self.conn.affected_rows(),
                           self.conn.insert_id())

    def selectElt(self, query):
        try:
            return self._query(query).getNextTuple()[0]
        except TypeError:
            raise Exception, "Query returned no records"

    def selectRec(self, query, fmt=MySQLResult.FMT_TUPLE):
        rec = self._query(query).getNextRec(fmt)
        if rec is None:
            raise Exception, "Query returned no records"
        else:
            return rec
    def selectTuple(self, query):
        return self.selectRec(query, fmt=MySQLResult.FMT_TUPLE)
    def selectDict(self, query):
        return self.selectRec(query, fmt=MySQLResult.FMT_DICT)
    def selectList(self, query):
        return self.selectRec(query, fmt=MySQLResult.FMT_LIST)
    def selectObject(self, query):
        return self.selectRec(query, fmt=MySQLResult.FMT_OBJECT)

    def selectAllRecs(self, query, fmt=MySQLResult.FMT_TUPLE):
        return self._query(query).getAllRecs(fmt)
    def selectAllTuples(self, query):
        return self.selectAllRecs(query, fmt=MySQLResult.FMT_TUPLE)
    def selectAllDicts(self, query):
        return self.selectAllRecs(query, fmt=MySQLResult.FMT_DICT)
    def selectAllLists(self, query):
        return self.selectAllRecs(query, fmt=MySQLResult.FMT_LIST)
    def selectAllObjects(self, query):
        return self.selectAllRecs(query, fmt=MySQLResult.FMT_OBJECT)

    def listDBs(self, like=None):
        query = 'SHOW DATABASES'
        if like is not None:
            query += ' LIKE %s' % self.stringLiteral(like)
        return self.selectField(query)
    def listTables(self, db=None, like=None):
        query = 'SHOW TABLES'
        if db is not None:
            query += ' FROM %s' % db
        if like is not None:
            query += ' LIKE %s' % self.stringLiteral(like)
        return self.selectField(query)
    def listFields(self, table, db=None, like=None):
        query = 'SHOW FIELDS FROM '
        if db is not None:
            query += '%s.' % db
        query += '%s' % table
        if like is not None:
            query += ' LIKE %s' % self.stringLiteral(like)
        return self.selectField(query)

    def selectField(self, query):
        recs = self.selectAllTuples(query)
        return map(operator.getitem, recs, [0] * len(recs))
    def selectAllFields(self, query, how=0):
        result = self._query(query)
        recs = result.getAllTuples()
        if len(recs) == 0: return None
        if how == 0:
            return zip(*recs)
        else:
            return dict(zip(result.getAllFieldNames(), zip(*recs)))

    def selectFieldsDict(self, query):
        recs = self.selectAllFields(query)
        d = {}
        map(operator.setitem, [d] * len(recs[0]), recs[0], recs[1])
        return d

    def buildSelect(self,
                    expr,
                    from_tables=None,
                    where=None,
                    group_by=None,
                    having=None,
                    order_by=None,
                    limit=None,
                    options=()
                    ):
        query = 'SELECT %s %s' % (" ".join(options), expr)
        if from_tables is not None:
            query += ' FROM %s' % from_tables
            if where is not None:
                query += ' WHERE %s' % where
            if group_by is not None:
                query += ' GROUP BY %s' % group_by
            if order_by is not None:
                query += ' ORDER BY %s' % order_by
            if limit is not None:
                if operator.isSequenceType(limit):
                    query += ' LIMIT %d, %d' % limit
                else:
                    query += ' LIMIT %d' % limit
        return query

    def selectCount(self, table, expr=None, distinct=0, where=None):
        if expr is None:
            expr = "*"
        elif distinct:
            expr = "DISTINCT %s" % expr

        return self.selectElt(self.buildSelect(expr=("COUNT(%s)" % expr),
                                               from_tables=table,
                                               where=where))

    def selectMedian(self, table, field, where=None):
        # Get total number of recs
        num_recs = self.selectCount(table, where)
        # If num recs is even, average of middle recs, else middle rec
        if num_recs % 2 == 0: # even, r
            vals = self.selectField(self.buildSelect(expr=field,
                                                     from_tables=table,
                                                     where=where,
                                                     limit=(num_recs / 2 - 1,
                                                            2)))
            return (vals[0] + vals[1]) / 2
        else: # odd
            return self.selectElt(self.buildSelect(expr=field,
                                                   from_tables=table,
                                                   where=where,
                                                   limit=num_recs / 2))
    def _selectFunc(self, func, table, field, where=None):
        return self.selectElt(self.buildSelect(expr='%s(%s)' % (func, field),
                                               from_tables=table,
                                               where=where))
    def selectAvg(self, table, field, where=None):
        return self._selectFunc('AVG', table, field, where)
    def selectMax(self, table, field, where=None):
        return self._selectFunc('MAX', table, field, where)
    def selectMin(self, table, field, where=None):
        return self._selectFunc('MIN', table, field, where)
    def selectStd(self, table, field, where=None):
        return self._selectFunc('STD', table, field, where)
    def selectSum(self, table, field, where=None):
        return self._selectFunc('SUM', table, field, where)
    def selectBitOr(self, table, field, where=None):
        return self._selectFunc('BIT_OR', table, field, where)
    def selectBitAnd(self, table, field, where=None):
        return self._selectFunc('BIT_AND', table, field, where)

    def update(self, table, data, where=None, limit=None):
        # If no data with which to update, just return
        if len(data) == 0: return
        query = 'UPDATE %s SET ' % table
        query += ', '.join(map(' = '.join,
                               zip(data.keys(),
                                   map(self.stringLiteral, data.values()))))
        if where is not None: query += ' WHERE %s' % where
        if limit is not None: query += ' LIMIT %d' % limit
        self._query(query)
        return self.conn.affected_rows()

    def _insertOrReplaceFields(self, cmd, table, fields, delayed=0):
        # if dictionary given, add field names clause
        if operator.isMappingType(fields):
            values = fields.values()
            fields = fields.keys()
        else:
            values = fields
            fields = None
        
        # determine number of recs to insert, if zero, just return
        numRecsList = map(len, filter(isTupleOrList, values))
        if len(numRecsList) == 0 or min(numRecsList) == 0:
            numRecs = 1
        else:
            numRecs = min(numRecsList)

        # adjust values that are not sequences to be sequences of
        # num_recs length containing the non sequence value
        # also clean values for the database
        filledVals = [None] * len(values)
        for i in range(len(values)):
            if not isTupleOrList(values[i]):
                filledVals[i] = ((values[i]),) * numRecs
            else:
                filledVals[i] = values[i]
        
        return self._insertOrReplaceRecs(cmd, table, zip(*filledVals),
                                         fields, delayed)

    def _insertOrReplaceRecs(self, cmd, table, recs,
                             fields=None, delayed=0):
        if delayed:
            query = '%s DELAYED INTO %s ' % (cmd, table)
        else:
            query = '%s INTO %s ' % (cmd, table)

        if len(recs) == 0:
            return

        # If we have a mapping type, extract the field
        # names from the first record
        if isinstance(recs[0], dict):
            fields = recs[0].keys()
            values = [[self.stringLiteral(rec[field])
                       for field in fields] for rec in recs]
        else:
            values = [[self.stringLiteral(elt)
                       for elt in rec] for rec in recs]
        if fields is not None:
            query += '(%s) ' % ', '.join(fields)
        query += 'VALUES (%s)' % '), ('.join(map(', '.join, values))
        self._query(query)
        return (self.conn.affected_rows(), self.conn.insert_id())
    def insertRec(self, table, rec, fields=None, delayed=0):
        self._insertOrReplaceRecs('INSERT', table, [rec], fields, delayed)
    def replaceRec(self, table, rec, fields=None, delayed=0):
        self._insertOrReplaceRecs('REPLACE', table, [rec], fields, delayed)
    def insertRecs(self, table, recs, fields=None, delayed=0): 
        self._insertOrReplaceRecs('INSERT', table, recs, fields, delayed)
    def replaceRecs(self, table, recs, fields=None, delayed=0):
        self._insertOrReplaceRecs('REPLACE', table, recs, fields, delayed)
    def insertFields(self, table, fields, delayed=0): 
        self._insertOrReplaceFields('INSERT', table, fields, delayed)
    def replaceFields(self, table, fields, delayed=0):
        self._insertOrReplaceFields('REPLACE', table, fields, delayed)
       
    def delete(self, table, where=None, order_by=None, limit=None):
        query = 'DELETE FROM %s' % table
        if where is not None: query += ' WHERE %s' % where
        if order_by is not None: query += ' ORDER BY %s' % order_by
        if limit is not None: query += ' LIMIT %d' % limit
        self._query(query)
        return self.conn.affected_rows()

    def truncate(self, table):
        self._query('TRUNCATE TABLE %s' % table)

    def dropTable(self, table, if_exists=1):
        query = 'DROP TABLE '
        if if_exists: query += 'IF EXISTS '
        if isTupleOrList(table):
            query += ', '.join(table)
        else:
            query += table
        self._query(query)

def test():
    db = MySQLDB(db='test')

    tableName = 'testDBMySQL'    
    db.dropTable(tableName)

    db.query("""
    CREATE TABLE %s (
        field01 TINYINT,
        field02 SMALLINT,
        field03 MEDIUMINT,
        field04 INT,
        field05 BIGINT,
        field06 FLOAT,
        field07 DOUBLE,
        field08 DECIMAL(10,2),
        field09 DATE,
        field10 DATETIME,
        field11 TIMESTAMP,
        field12 TIME,
        field13 YEAR,
        field14 CHAR(10),
        field15 VARCHAR(10),
        field16 TINYTEXT,
        field17 TEXT,
        field18 MEDIUMTEXT,
        field19 LONGTEXT,
        field20 ENUM('value1', 'value2', 'value3'),
        field21 SET('value1', 'value2', 'value3')
    )""" % tableName)

    from time import strftime

    data = {'field01': (1, 2, 3, 4),
            'field02': (10, 20, 30, 40),
            'field03': (100, 200, 300, 400),
            'field04': (1000, 2000, 3000, 4000),
            'field05': (10000, 20000, 30000, 40000),
            'field06': (1.0, 2.0, 3.0, 4.0),
            'field07': (10.0, 20.0, 30.0, 40.0),
            'field08': (1.23, 4.56, 7.89, 10.23),
            'field09': (strftime('%Y-%m-%d'),) * 4,
            'field10': (strftime('%Y-%m-%d %H:%M:%S'),) * 4,
            'field11': (strftime('%Y-%m-%d %H:%M:%S'),) * 4,
            'field12': (strftime('%H:%M:%S'),) * 4,
            'field13': (strftime('%Y'),) * 4,
            'field14': ('value1', 'value2', 'value3', 'value4'),
            'field15': ('value1', 'value2', 'value3', 'value4'),
            'field16': ('value1', 'value2', 'value3', 'value4'),
            'field17': ('value1', 'value2', 'value3', 'value4'),
            'field18': ('value1', 'value2', 'value3', 'value4'),
            'field19': ('value1', 'value2', 'value3', 'value4'),
            'field20': ('value1', 'value2', 'value3', 'value4'),
            'field21': ('value1', 'value2', 'value3', 'value4'),
            }
    db.insertFields(tableName, data)

    query = "SELECT * FROM %s" % tableName
    
    result = db.query(query)
    print result.getFieldName(0)
    print result.getFieldType(0)
    print result.getFieldLen(0)
    print result.getAllFieldNames()
    print result.getAllFieldTypes()
    print result.getAllFieldLengths()
    print result.getNextTuple()
    print result.getNextDict()
    print result.getNextList()
    print result.getNextObject()
    result.seekToRec(0)
    print result.getAllTuples()
    result.seekToRec(0)
    print result.getAllDicts()
    print result.getNumRecs()
    print result.getNumFields()

    print result.getQuery()
    print result.getNumAffectedRecs()
    print result.getInsertID()

    print db.listDBs()
    print db.listTables()
    print db.listFields(tableName)

    print db.selectElt(query)
    print db.selectTuple(query)
    print db.selectDict(query)
    print db.selectList(query)
    print db.selectObject(query)
    print db.selectAllRecs(query)
    print db.selectAllTuples(query)
    print db.selectAllDicts(query)
    print db.selectAllLists(query)
    print db.selectAllObjects(query)
    print db.selectField(query)
    print db.selectAllFields(query)
    print db.selectAllFields(query, how=1)
    print db.selectFieldsDict(query)

    db.insertRecs(tableName, recs=db.selectAllObjects(query))
    print db.selectCount(tableName)
    db.insertRecs(tableName, recs=db.selectAllDicts(query))
    print db.selectCount(tableName)
    db.insertRecs(tableName, recs=db.selectAllTuples(query))
    print db.selectCount(tableName)

    for rec in db.selectAllFields(query):
        print rec
    
    db.dropTable(tableName)

def time(num):
    db = MySQLDB(db='bandroster')
    from time import time
    start = time()
    for i in range(num):
        db.listFields('user', db='bandroster', like='u%')
    end = time()
    return end - start
