from argparse import Namespace
import pytest
from bifrostlib import common
from bifrostlib import datahandling
from bifrostlib import database_interface
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import RunReference
from bifrostlib.datahandling import Run
from bifrost_salmonella_subspecies_dtartrate import launcher
import pymongo
import os
import shutil
from pathlib import Path

@pytest.fixture
def test_connection():
    assert datahandling.has_a_database_connection()
    assert "TEST" in os.environ['BIFROST_DB_KEY'].upper()  # A very basic piece of protection ensuring the word test is in the DB

class TestBifrostSubspeciesDtartrate:
    component_name = "salmonella_subspecies_dtartrate__v1.0.5"
    bifrost_install_dir = Path(os.environ['BIFROST_INSTALL_DIR'])
    bifrost_config_and_data_path = Path(f"{bifrost_install_dir}/bifrost/test_data")
    current_dir = os.getcwd()
    test_dir = bifrost_config_and_data_path / "output/test__salmonella_subspecies_dtartrate"
    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"}, 
            "name": "SRR2094561", 
            "components": [], 
            "categories": {
                "paired_reads": {
                    "summary": {
                        "data": [str(bifrost_config_and_data_path/"/samples/SRR2094561_1.fastq.gz"),
                                 str(bifrost_config_and_data_path/"/samples/SRR2094561_2.fastq.gz")]
                    }
                },
                "mlst": {
                    "summary": {
                        "sequence_type": {
                            "senterica":"34"
                        }
                    }
                },
                "species_detection": {
                    "summary": {
                        "species": "Salmonella enterica"
                    }
                }
            }
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]

    @classmethod
    def setup_class(cls):
        client = pymongo.MongoClient(os.environ['BIFROST_DB_KEY'])
        db = client.get_database()
        cls.clear_all_collections(db)
        col = db["samples"]
        col.insert_many(cls.bson_entries)
        launcher.initialize()
        os.chdir(cls.current_dir)

    @classmethod
    def teardown_class(cls):
        client = pymongo.MongoClient(os.environ['BIFROST_DB_KEY'])
        db = client.get_database()
        cls.clear_all_collections(db)

    @staticmethod
    def clear_all_collections(db):
        db.drop_collection("components")
        db.drop_collection("hosts")
        db.drop_collection("run_components")
        db.drop_collection("runs")
        db.drop_collection("sample_components")
        db.drop_collection("samples")

    def test_info(self):
        launcher.main(["--info"])

    # def test_help(self):
    #     launcher.main(["--help"])

    def test_pipeline(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

        os.mkdir(self.test_dir)
        test_args = [
            "--sample_name", "SRR2094561",
            "--outdir", str(self.test_dir)
        ]
        launcher.main(args=test_args)
        assert os.path.isfile(f"{self.test_dir}/{self.component_name}/datadump_complete")
        #shutil.rmtree(self.test_dir)
        #assert not os.path.isdir(f"{self.test_dir}/{self.component_name}")


