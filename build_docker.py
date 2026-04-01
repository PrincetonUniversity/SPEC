#!/usr/bin/env python3
import docker
import os
from datetime import datetime

###### THIS IS JUST FOR PUSHING TO SHARED REPO
dockerfile = 'Dockerfile'
image_name = 'spec'
image_tag = datetime.now().strftime("%m%d%y")
platform   = 'linux/amd64'
remote_repository = 'containers.qarnot.com'
remote_user = 'qrntadgaussf'
remote_username = 'qrnt_adm@gauss-fusion.com'
build_args = ['platform',platform]
remote_password = os.getenv('QARNOT_ADM_REGISTRY_KEY')



# Create the Client
client = docker.from_env()
# Login to the remote Server
client.login(username=remote_username,registry=remote_repository,password=remote_password)
# Build the image
(image,log)=client.images.build(path='./',tag=f"{image_name}:{image_tag}",platform=platform,dockerfile=dockerfile)
# Now tag the image with the remote name and 'latest'
remote_name = f'{remote_repository}/{remote_user}/{image_name}'
late_tag = 'latest'
image.tag(remote_name,image_tag)
image.tag(remote_name,late_tag)
# Push the image tags
for line in client.images.push(remote_name, tag=image_tag,stream=True,decode=True):
	print(line)
for line in client.images.push(remote_name, tag=late_tag,stream=True,decode=True):
	print(line)
