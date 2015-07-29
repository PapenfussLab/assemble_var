"""
SQLAlchemy tables for mungo classes
"""

from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker, mapper
from sqlalchemy import Table, Column, Integer, Float, Text, String


class Alchemy(object):
    def __init__(self, engine=None):
        self.engine = None
        self.Session = None
        self.metadata = MetaData()
    
    def startSession(self, dsn):
        self.engine = create_engine(dsn)
        self.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine, autoflush=True, transactional=True)
        return self.Session()
