"""This module helps connect to a PostgreSQL database"""
from dotenv import load_dotenv
import psycopg2
import sys
import os

# Use the parameters in database.ini to configure the database connection
def config():
  """Collect credentials for connecting to database using environment variables 

    This function parsers out the PostgreSQL credentials from environment

    :return: database connection credentials
    :rtype: dict
    :raises: :exc:`FileNotFound`

    :example .env:
      .. code-block:: env
        database=database_name
        user=database_owner
        password=password
        host=localhost
        port=5432
  """
  load_dotenv()
  # get section, default to postgresql
  db = {}
  required_sections = [ 'database', 'user', 'password', 'host', 'port' ]
  for section in required_sections:
    try:
      db[section] = os.getenv(section)
    except:
      raise Exception(f"Section '{section}' not found. Please define it as an environment variable.")

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
