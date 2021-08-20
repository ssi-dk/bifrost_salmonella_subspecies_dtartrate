# bifrost_salmonella_subspecies_dtartrate

This component is run given a sample id already added into the bifrostDB. If the sample is registered as Salmonella enterica, this will look up the serotype in enterobase using the ST.

## Howto launch
  /bifrost/components/bifrost_salmonella_subspecies_dtartrate$ docker-compose run bifrost_salmonella_subspecies_dtartrate
## How to debug
  /bifrost/components/bifrost_salmonella_subspecies_dtartrate$ docker-compose run --entrypoint=bash bifrost_salmonella_subspecies_dtartrate
