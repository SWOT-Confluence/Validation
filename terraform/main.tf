 terraform {
  backend "s3" {
    encrypt = true
  }
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  default_tags {
    tags = local.default_tags
  }
  region  = var.aws_region
}

data "aws_caller_identity" "current" {}

locals {
  account_id = sensitive(data.aws_caller_identity.current.account_id)
  default_tags = length(var.default_tags) == 0 ? {
    application : var.app_name,
    environment : lower(var.environment),
    version : var.app_version
  } : var.default_tags
}

module "confluence-validation" {
  source            = "./modules/validation"
  app_name          = var.app_name
  app_version       = var.app_version
  aws_region        = var.aws_region
  environment       = var.environment
  prefix            = var.prefix
}