from sqlalchemy import create_engine, MetaData
from sqlalchemy import Table, Column, Integer, Text, Float
from sqlalchemy.orm import sessionmaker


alchemyTypes = {
    int: Integer,
    float: Float
}


def startSession(dsn):
    """Create a database (if necessary) and return a session object
    
    @param dsn: database source name (DSN)
    @param metadata: SQLAlchemy metadata
    @returns: session
    """
    meta = Meta()
    meta.engine = create_engine(dsn)
    meta.metadata.create_all(engine)
    meta.Session = sessionmaker(autoflush=True, transactional=True)
    meta.Session.configure(bind=engine)
    session = meta.Session()
    return session


def createTable(tableName, metadata, attributes, converterList, indexedAttributes=None):
    """Create table from attribute and converter lists
    
    @param tableName: database table name
    @param metadata: SQLAlchemy metadata
    @param attributes: attribute list
    @param converterList: converter list; [(attribute, type), ...]; 
    default type is str
    @returns: table
    """
    converters = dict(converterList)
    if not indexedAttributes: indexedAttributes = []
    table = Table(tableName, metadata)
    for attrib in attributes:
        if converters.has_key(attrib):
            table.append_column(
                Column(attrib, alchemyTypes[converters[attrib]], 
                index=attrib in indexedAttributes))
        else:
            table.append_column(Column(attrib, Text, 
                index=attrib in indexedAttributes))
    table.append_column(Column('id', Integer, primary_key=True))
    return table

