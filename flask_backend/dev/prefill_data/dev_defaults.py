# config/test_config.py
from flask_backend.dev.tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ, TEST_SEQ_1, TEST_SEQ_2

dev_config = {
    "templateSequence": TEST_TEMPLATE_SEQ,
    "species": "",
    "kozak": "",
    "maxMutationsPerSite": "one",
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

dev_config_2 = {
    "templateSequence": TEST_TEMPLATE_SEQ,
    "species": "",
    "kozak": "",
    "maxMutationsPerSite": "one",
    "verbose_mode": True,
    "sequencesToDomesticate": [
        {
            "sequence": TEST_SEQ_1,
            "primerName": "NLS-mTag-BFP",
            "mtkPartLeft": "5",
            "mtkPartRight": "5"
        },{
            "sequence": TEST_SEQ_2,
            "primerName": "hygroR",
            "mtkPartLeft": "6",
            "mtkPartRight": "6"
        },
    ]
}