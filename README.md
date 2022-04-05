# Docker 101

## Pre-requisites

Required reading:

- https://staphb.org/docker-builds/run_containers/
- https://staphb.org/docker-builds/make_containers/

## Hello world

> This basic container uses the Alpine Linux distribution and without any modification simply echoes "Hello world" to the host machine. Objectives:
> - Building a container
> - Running a container with a default `CMD`
> - Overriding the default `CMD`
> - Stepping inside the container with `-it ... /bin/ash`

Build:

```sh
docker build -t d101/helloworld ./01-helloworld
```

Run:

```sh
docker run --rm d101/helloworld
docker run --rm d101/helloworld echo "Overriding the built in cmd with something else"
docker run --rm d101/helloworld whoami
docker run --rm d101/helloworld pwd
docker run --rm d101/helloworld date
docker run --rm d101/helloworld hostname
```

Step inside the running container

```sh
docker run --rm -it d101/helloworld /bin/ash
whoami
pwd
ls -l
hostname
exit
```

## HTSlib

> Build a container installing several required dependencies, then download and install bcftools. Objectives:
> - Build a container with external software and dependencies
> - Run the container
> - Mount host volumes to the running container with `-v`
> - Run tools inside the container on data on the host system via a volume mount
> - Step inside the running container and play around

Build:

```sh
docker build -t d101/htslib ./02-htslib
```

Run:

```sh
docker run --rm d101/htslib
docker run --rm d101/htslib samtools
docker run --rm d101/htslib wgsim
docker run --rm d101/htslib bcftools
docker run --rm -v $(pwd)/exampledata:/data d101/htslib bcftools view /data/a.vcf.gz
```

Hop inside the running container with the example data mounted to `/data/`:

```sh
# Start and step inside container
docker run -it --rm -v $(pwd)/exampledata:/data d101/htslib /bin/ash

# Inside the container!
ls -l
whoami
pwd
bcftools merge a.vcf.gz b.vcf.gz --force-samples | bcftools query -f "%CHROM\t%POS\t[%TGT\t]\n"
exit
```

## Bcftools as an executable

> This demonstrates a different way to use Docker, essentially using it as an executable with an `ENTRYPOINT` in the Dockerfile. Objectives:
> - Mount the current host pwd to the same pwd in the container
> - Set the working dir in the container to the same full path on the host
> - Run a container as an executable on host data

Build:

```sh
docker build -t d101/bcftools ./03-bcftools
```

Run:

```sh
cd exampledata
docker run --rm d101/bcftools
docker run --rm -v $(pwd):$(pwd) -w $(pwd) d101/bcftools view b.vcf.gz
docker run --rm -v $(pwd):$(pwd) -w $(pwd) d101/bcftools query -f '[%GT]\n' b.vcf.gz | sort | uniq -c
docker run --rm -v $(pwd):$(pwd) -w $(pwd) d101/bcftools query -f '[%GT]\n' b.vcf.gz | sort | uniq -c
cd ..
```

## Pipeline

> Demonstrates how to build an entire pipeline that runs inside a container using data from the host system inside the container. Objectives:
> - Set the user and group ID of the user _inside the container_
> - Mount a host directory to a mount point on the container
> - Pass an environment variable to the running container
> - Show how the container runs several scripts on data on the host

Build:

```sh
docker build -t d101/pipeline ./04-pipeline
```

Run:

```sh
vcf_dir=$(pwd)/exampledata

docker run -u $(id -u):$(id -g) --rm -v $vcf_dir:/data -e SAMPLE="a" d101/pipeline

ls -lh exampledata
cat exampledata/a.out.csv
open exampledata/a.out.csv
```

## Learn more

- Inspect Dockerfiles at https://github.com/stephenturner/slimbioinfo
- Inspect Dockerfiles at https://github.com/StaPH-B/docker-builds
