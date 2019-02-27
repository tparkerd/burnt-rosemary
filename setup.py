if __name__=="__main__":
  from setuptools import setup
  setup(
    name="PrelimInsertionPipeline",
    version="1.0-alpha",
    install_requires=[
      "pandas",
      "psycopg2-binary",
      "configparser",
      "tqdm"
    ]
  )
