# config/test_config.py
from flask_backend.tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
CONFIG = {
    "templateSequence": TEST_TEMPLATE_SEQ,
    "species": "",
    "kozak": "",
    "max_mut_per_site": 1,
    "verbose_mode": True,
    "sequencesToDomesticate": [
        {
            "sequence": TEST_SEQ,
            "primerName": "NLS-mTag-BFP-hygroR",
            "mtkPartLeft": "6",
            "mtkPartRight": "6"
        }
    ]
}
