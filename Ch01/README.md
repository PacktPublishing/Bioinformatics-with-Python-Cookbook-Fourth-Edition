This is the README for Ch01.


These are the commands to build and run the Docker container for the book:

docker build -t bio https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-fourth-edition.git#main:docker/main

docker run -ti -p 9875:9875 -v /Users/shanebrubaker/work/docker_files:/data bio 


