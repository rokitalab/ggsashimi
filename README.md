# ggsashimi

Command-line tool for the visualization of splicing events across multiple samples, adapted for use with rokita-lab and mounting cavatica project bams on an EC2 instance. Utilizes a modified `ggsashimi.py` python script from guigolab/ggsashimi for plot generation. 

## Installation<a name="installation"></a>

From an EC2 instance:

1. Clone the repository:
```
git clone git@github.com:rokitalab/ggsashimi.git
```

2. Pull Docker container:
```
docker pull pgc-images.sbgenomics.com/rokita-lab/ggsashimi:latest
```

3. Start the Docker container from the `ggsashimi` folder:
```
docker run --privileged --name <NAME> -d -e PASSWORD=<ANYTHING> -p 80:8787 -v $PWD:/home/rstudio/ggsashimi pgc-images.sbgenomics.com/rokita-lab/ggsashimi:latest
```
To launch Rstudio in your browser, navigate to the instance IP address in your web browser. The username for login is rstudio and the password is set in the docker run command above.

Or access the docker container via command line 
```
docker exec --privileged -ti <NAME> bash
```

4. Configure SBFS credentials: 
Run `sbfs configure` and enter the following when prompted:
API endpoint [None]: `https://cavatica-api.sbgenomics.com/v2`
Authentication token [None]: (personal CAVATICA authentication token)
NOTE: these parameters will automatically be assigned as the “default” profile in the configuration file.

5. Download reference files (GRCh38 genome fasta and GENCODE v39 GTF)

```
bash download_data.sh
```

6. Run ggsashimi

Generate sashimi plots; for example:

```
bash run_ggsashimi.sh --sample_file examples/samples.txt --coord_file examples/regions.txt
```

## How to run ggsashimi

`run_ggsashimi.sh` requires two arguments: a `sample_file` and `coord_file`

* `sample_file` must contain the following columns: 

1. `sample_name`: name to be used in bam mapping file
2. `cavatica_project`: project to mount and pull cram files
3. `cram_name`: name of sample cram file to be pulled from cavatica project
4. `group`: sashimi plot group label

see `examples/samples.txt` for formatting:

| sample_name          | cavatica_project               | cram_name                                                         | group               |
|----------------------|--------------------------------|-------------------------------------------------------------------|---------------------|
| GTEx-cerebellum-1    | sicklera/pbta-and-normal-crams | b7bccab5-2f6f-48bd-b9c1-54a27be86ebd.Aligned.out.sorted.cram      | GTEx - Cerebellum   |
| GTEx-cerebellum-2    | sicklera/pbta-and-normal-crams | 041b942c-4fc7-4a3c-aa56-64c53b0da2d9.Aligned.out.sorted.cram      | GTEx - Cerebellum   |
| GTEx-skin-1          | sicklera/pbta-and-normal-crams | 9e00d670-eb1a-4d43-8f8c-2eb10cba020a.Aligned.out.sorted.cram      | GTEx - Skin         |
| GTEx-skin-2          | sicklera/pbta-and-normal-crams | d83776bb-ac69-43d0-825d-e56db280888d.Aligned.out.sorted.cram      | GTEx - Skin         |
| GTEx-testis-1        | sicklera/pbta-and-normal-crams | aa473c90-cd79-41b2-a72d-d4b4de27d7b8.Aligned.out.sorted.cram      | GTEx - Testis       |
| GTEx-testis-2        | sicklera/pbta-and-normal-crams | 8297f79b-ed25-4897-9201-3e67d732fb0c.Aligned.out.sorted.cram      | GTEx - Testis       |


* `coord_file` must contain the following columns:

1. `region`: coordinates in chr:start-end format to plot. It is recommended to include at least 500bp on either side of interval of interest for sashimi plot generation
2. `name`: name for tmp file prefixes and output plot file prefix. 

see `examples/regions.txt` for formatting: 

| region                    | name  |
|---------------------------|-------|
| chr2:200859569-200861516  | CLK1  |
| chr1:151324500-151327600  | PI4KB |
| chr7:108234500-108240250  | NRCAM |
| chr22:20993000-20997000   | LZTR1 |

Other optional arguments:

* `--min_coverage`: Minimum number of reads supporting a junction to be drawn [default=10]

