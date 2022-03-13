## Usage

- Change parameters in *run.sh* according to your needs

```shell
cd folder
./run.sh
```

## Devcloud
Submit job
```shell
$ ./q devcloud-run
```

[VTune on devcloud](VTune
https://www.intel.com/content/www/us/en/develop/documentation/vtune-cookbook/top/configuration-recipes/using-vtune-server-with-vs-code-intel-devcloud.html)
Advisor on devcloud 
```shell
$ qsub -I
$ advisor -collect=roofline -project-dir=./roofline <your-programm> <programm-args>
$ advixe-cl -report=roofline -project-dir=./roofline -report-output ./roofline.html
```
Then you can download *roofline.html* using _scp_ or other tools