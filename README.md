# validation

Validation serves as a way to validate results in the Confluence workflow. It compares gage data and produces stats and hydrographs for results that match both in reach identifier and temporally.

## installation

## setup

## execution

## deployment

There is a script to deploy the Docker container image and Terraform AWS infrastructure found in the `deploy` directory.

Script to deploy Terraform and Docker image AWS infrastructure

REQUIRES:

- jq (<https://jqlang.github.io/jq/>)
- docker (<https://docs.docker.com/desktop/>) > version Docker 1.5
- AWS CLI (<https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html>)
- Terraform (<https://developer.hashicorp.com/terraform/tutorials/aws-get-started/install-cli>)

Command line arguments:

[1] registry: Registry URI
[2] repository: Name of repository to create
[3] prefix: Prefix to use for AWS resources associated with environment deploying to
[4] s3_state_bucket: Name of the S3 bucket to store Terraform state in (no need for s3:// prefix)
[5] profile: Name of profile used to authenticate AWS CLI commands

Example usage: ``./deploy.sh "account-id.dkr.ecr.region.amazonaws.com" "container-image-name" "prefix-for-environment" "s3-state-bucket-name" "confluence-named-profile"`

Note: Run the script from the deploy directory.
