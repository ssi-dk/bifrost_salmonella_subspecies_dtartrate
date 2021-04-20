from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import Category
from typing import Dict
import os

def extract_serotype_results(serotype: Category, results: Dict, component_name: str) -> None:
    file_name = "serotype.txt"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)

    for line in open(file_path,'r'):
        (ID, ST, serotype1, count1, serotype2, count2) = line.strip().split("\t")
    results['Enterobase_serotype1'] = serotype1
    results['Enterobase_count1'] = count1
    results['Enterobase_serotype2'] = serotype2
    results['Enterobase_count2'] = count2

    serotype["summary"]["serotype"] = results["Enterobase_serotype1"]


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
                "serotype": "",
                "antigenic profile": "",
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