#!/usr/bin/env python3

import boto3

from ec2_metadata import ec2_metadata
instance_id = ec2_metadata.instance_id
ec2 = boto3.resource('ec2')
instance = [i for i in ec2.instances.all() if i.id == instance_id][0]
name_tag = [tag for tag in instance.tags if tag['Key'] == 'Name'][0]
return name_tag['Value']
