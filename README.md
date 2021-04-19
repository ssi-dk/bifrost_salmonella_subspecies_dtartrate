# bifrost_enterobase

This component is run given a sample id already added into the bifrostDB. If the sample is registered as Salmonella enterica, this will look up the serotype in enterobase using the ST.

## Howto launch
  /bifrost/components/bifrost_enterobase$ docker-compose run bifrost_enterobase
## How to debug
  /bifrost/components/bifrost_enterobase$ docker-compose run --entrypoint=bash bifrost_enterobase
