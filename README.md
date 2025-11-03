# ggsashimi

Command-line tool for the visualization of splicing events across multiple samples, adapted for use with rokita-lab and mounting cavatica project bams on an EC2 instance.

## Installation<a name="installation"></a>

From an EC2 instance:

1. Clone the repository:
```
git clone git@github.com:rokitalab/splicing-variants.git
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

5. Mount a cavatica project to the data directory:
```
mkdir data/cavativa
sbfs mount --profile default --project sicklera/pbta-and-normal-crams data/cavatica
```

To unmount, run `sbfs unmount <mount_dir>`
