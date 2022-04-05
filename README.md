# docker101

## Hello world

Build:

```sh
docker build --no-cache -t d101/helloworld ./01-helloworld
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

## HTSlib

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