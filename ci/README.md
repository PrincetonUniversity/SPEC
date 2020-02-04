# Simple gitlab CI notes

## Docker Containers

Docker is a tool which helps us control our envirnment with code.
A "container" is like a virtual machine, but lighter,
because it shares some features with the host system.
They also have a hashing scheme which allows re-use of filesystem "layers",
if docker beleives it is safe to do so.  So for example,
if you have 3 containers that you use for different projects and they all
start with centos7 base, that will probably be cached only once.

I will assume you have Docker installed and working for linux development.
It is pretty straightforward these days.

### Dockerfile

This is a file which sequentially builds up a container.  We usually
start from a known base, and add the minimal set of things we need.
It is good practice to control version, but I do not always.  There
are many excellent (also some poor, so watch out) guides online for
good practices.

Here we'll build up a Dockerfile for `SPEC`.

First we start with some base image:

```
FROM centos:7
```

Then we'll add things we need

```
#setup the base sysem
RUN yum group install "Development Tools" -y
RUN yum -y install epel-release
RUN yum install gcc-gfortran openmpi openmpi-devel hdf5 hdf5-devel -y
```

And from experience, we have more config to get openmpi working:

```
# Complete the setup of openmpi, they used environmentmodules, blehgh
ENV PATH="/usr/lib64/openmpi/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/lib64/openmpi/lib:${LD_LIBRARY_PATH}"
```

I save that all in a file with a simple name, like Dockerfile.centos7.

It is important for posterity that we save it with the version controlled code.

### Building a Docker container

```
docker build -t my_tag_name -f Dockerfile.centos7 .
```

Tags are simply a way to name/label/version things.
Below I use testspec.

The first build on any docker host will usually take some time.
Docker will need to download lots of things,
and also go through the package installation processes.

Future passes will only rebuild changesets. If you don't change anything,
usually docker will use its local cache.  Here you can see evidence of the
hashes for the layers.

```
# docker build -t testspec -f Dockerfile.centos7 .
Sending build context to Docker daemon  2.048kB
Step 1/6 : FROM centos:7
 ---> 9f38484d220f
Step 2/6 : RUN yum group install "Development Tools" -y
 ---> Using cache
 ---> e98fc3dd0ef4
Step 3/6 : RUN yum -y install epel-release
 ---> Using cache
 ---> 3da7096c3718
Step 4/6 : RUN yum install gcc-gfortran openmpi openmpi-devel hdf5 hdf5-devel -y
 ---> Using cache
 ---> 301dbe27c50e
Step 5/6 : ENV PATH="/usr/lib64/openmpi/bin:${PATH}"
 ---> Using cache
 ---> 6c5a9aa04026
Step 6/6 : ENV LD_LIBRARY_PATH="/usr/lib64/openmpi/lib:${LD_LIBRARY_PATH}"
 ---> Using cache
 ---> 5649fc4ec484
Successfully built 5649fc4ec484
Successfully tagged testspec:latest
```


### Running a container locally

Sometimes containers are setup to have the code/application installed already.
That would be a `deployment`.

For development, we will want to access our local version of the code.
This requires mounting/binding your working directory.

This is one way to do that, if you are in the directory with your code:

```
docker run -v `pwd`:`pwd` -w `pwd` my_tag_name yourcommand
```

Try it:

```
➜  testspec docker run -v `pwd`:`pwd` -w `pwd` testspec uname -a
Linux 77f967b2d9b8 4.19.76-linuxkit #1 SMP Thu Oct 17 19:31:58 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
➜  testspec uname -a
Darwin gwright-lt 18.7.0 Darwin Kernel Version 18.7.0: Thu Jan 23 06:52:12 PST 2020; root:xnu-4903.278.25~1/RELEASE_X86_64 x86_64
```

Neat right?

### Pushing

So usually containers are stored in a `registry`.
These can be privately managed or public.  Dockerhub is the most notable
one.  Gitlab.com also offers a docker container registry for open source
projects.  Here I will demonstrate pushing the container to my dockerhub.
To use a registry like dockerhub, you need an account (its free to make one).
You will need to be logged in.  There are ways to store your credentials, but
the command to use uname/passwd is:

```
docker login
```

First we tag it again, now adding the registry to the tagname. Then we push.
If this is the first time pushing large layers, it may take a long time to upload.

```
docker tag testspec garrettwrong/testspec
docker push garrettwrong/testspec
```

Cool, now any docker with internet access can pull this image from the registry.


## Gitlab CI

Gitlab CI is a stupid simple approach to CI.  It uses a `yaml` file to store
shell commands, and also has syntactic sugar for really common CI operations.
Most of the other CI tools have ripped this off by now,
(travis, circle, github actions).
Gitlab should consider it flattering, they are still by far the best...
but you can use whichever one makes sense for your organization.
Here we have in internal gitlab that I am trying to get integration enabled for...

Some things you can do pretty easily:

* Use containers to build/run against multiple controlled environments
* Use something called `gitlab runners` to run on existing machiens (assuming you have permissions etc)
* Store results for a set period of time, or perform differing tasks on failure.
* Run tasks in pipelines stages in parallel
* Get automatic emails/feedback

The basic ci script stubbed here would do most of that already,
once its filled in of course.

You may read the .gitlab-ci.yml in the root of the git directory.

## Putting it together

So now we have a tool to make our own environments with code.

We have a tool to run arbitrary commands with our specified environments.

If we add just a little bit of structure, we'll have a pretty simple and
effective end to end automation.  Again, this is all described in controlled
code that is mated with your scientific code...

A simple structure I have used for several similar physics projects is to
create a directory structure of `ci/case_names`. For each case I will include
input, configuration, and reference output files.

If these files are large, or will change often, it might be advisable to store in a seperate repo,
or git submodule repo.  This is because it can cause your repo to get bloated in size.
