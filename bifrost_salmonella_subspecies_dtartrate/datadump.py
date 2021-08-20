from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import Category
from typing import Dict
import os

def extract_serotype_results(serotype: Category, results: Dict, component_name: str) -> None:
    file_name = "subspecies_dtartrate.txt"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)

    for line in open(file_path,'r'):
        (ID, ST, dtartrate, subspecies) = line.strip().split("\t")
    results['Subspecies'] = subspecies
    results['D-tartrate_pos10'] = dtartrate

    serotype["summary"]["Subspecies"] = results["Subspecies"]
    serotype["summary"]["D-tartrate_pos10"] = results["D-tartrate_pos10"]


def datadump(samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    component = Component.load(samplecomponent.component)
    
    serotype = samplecomponent.get_category("serotype")
    if serotype is None:
        serotype = Category(value={
            "name": "serotype",
            "component": samplecomponent.component,
            "summary": {
                "Subspecies": "",
                "D-tartrate_pos10": "",
            },
            "report": {}
        })
    extract_serotype_results(serotype, samplecomponent["results"], samplecomponent["component"]["name"])
    samplecomponent.set_category(serotype)
    sample.set_category(serotype)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    with open(os.path.join(samplecomponent["component"]["name"], "datadump_complete"), "w+") as fh:
        fh.write("done")


datadump(
    snakemake.params.samplecomponent_ref_json,
)
