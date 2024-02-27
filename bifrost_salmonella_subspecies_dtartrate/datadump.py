import os
from typing import Dict
from pprint import pprint

from bifrostlib import common
from bifrostlib.datahandling import (Category, Component, Sample,
                                     SampleComponent, SampleComponentReference)


def extract_subspecies_results(serotype: Category, results: Dict, component_name: str, file_name: str) -> None:
    #file_path = os.path.join(component_name, file_name)
    for line in open(file_name,'r'):
        subspecies = line.strip().split("\t")
    results['Subspecies'] = subspecies
    serotype["summary"]["Subspecies"] = results["Subspecies"]

def extract_dtartrate_results(serotype: Category, results: Dict, component_name: str, file_name: str) -> None:
    for line in open(file_name,'r'):
        dtartrate = line.strip()
    results['D-tartrate_pos10'] = dtartrate
    serotype["summary"]["D-tartrate_pos10"] = results["D-tartrate_pos10"]


def datadump(input: object, output: object, samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    component = Component.load(samplecomponent.component)
    
    serotype = sample.get_category("serotype")
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
    extract_subspecies_results(
        serotype,
        samplecomponent["results"],
        samplecomponent["component"]["name"],
        input.subspecies)
    extract_dtartrate_results(
        serotype,
        samplecomponent["results"],
        samplecomponent["component"]["name"],
        input.dtartrate)
    samplecomponent.set_category(serotype)
    sample.set_category(serotype)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    pprint(output.complete)
    with open(output.complete[0], "w+") as fh:
        fh.write("done")


datadump(
    snakemake.input,
    snakemake.output,
    snakemake.params.samplecomponent_ref_json,
)
