"""
Cloud API clients for variant annotation.
Provides cloud-based access to VEP, ClinVar, Ensembl, and gnomAD APIs.
"""

from cloud_api.base_client import BaseCloudClient
from cloud_api.clinvar_client import ClinVarEutilsClient
from cloud_api.vep_client import VEPCloudClient
from cloud_api.ensembl_client import EnsemblCloudClient
from cloud_api.gnomad_client import GnomADConstraintClient

__all__ = [
    "BaseCloudClient",
    "ClinVarEutilsClient",
    "VEPCloudClient",
    "EnsemblCloudClient",
    "GnomADConstraintClient",
]
