"""This module helps connect to a PostgreSQL database"""
from configparser import ConfigParser
import psycopg2
import sys

# Use the parameters in database.ini to configure the database connection
def config(filename='database.ini', section='postgresql'):
  """Parses ini file 

    This function parsers out the PostgreSQL credentials from input .ini

    :param filename: input file
    :type filename: str
    :param section: section within .ini
    :type j: str
    :return: database connection credentials
    :rtype: dict
    :raises: :exc:`FileNotFound`

    :example .ini:
      .. code-block:: ini

        [postgresql]
        database=database_name
        user=database_owner
        password=password
        host=localhost
        port=5432
  """
  # create a parser
  parser = ConfigParser()
  # read config file
  parser.read(filename)

  # get section, default to postgresql
  db = {}
  if parser.has_section(section):
    params = parser.items(section)
    for param in params:
      db[param[0]] = param[1]
  else:
    raise Exception('Section {0} not found in the {1} file'.format(section, filename))

  return db

# Return a connection to the database
def connect():
  """Creates connection object to database 

    This function creates a connection object to database

    :return: connection to database
    :rtype: connection object
  """

  conn = None
  try:
    # read connection parameters
    params = config()

    # connect to the PostgreSQL server
    print('Connection to the PostgreSQL database...')
    conn = psycopg2.connect(**params)

    # create a cursor
    cur = conn.cursor()

    print('PostgreSQL database version:')
    cur.execute('SELECT version()')

    db_version = cur.fetchone()
    print(db_version)

    # Close the connection
    cur.close()
  except (Exception, psycopg2.DatabaseError) as error:
    print('Unable to connect!\n{0}').format(error)
    sys.exit(1)
  
  return conn
