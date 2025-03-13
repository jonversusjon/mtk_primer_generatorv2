# config/test_config.py

from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
CONFIG = {
    "development": {
        "templateSequence": TEST_TEMPLATE_SEQ,
        "species": "",
        "kozak": "",
        "max_mut_per_site": 1,
        "verbose_mode": True,
        "sequencesToDomesticate": [
            {
                "sequence": TEST_SEQ,
                "primerName": "Test_1",
                "mtkPartLeft": "6",
                "mtkPartRight": "6"
            },
            {
                "sequence": "",
                "primerName": "",
                "mtkPartLeft": "",
                "mtkPartRight": ""
            }
        ]
    },
    "testing": {
        "templateSequence": "",
        "species": "",
        "kozak": "",
        "max_mut_per_site": 1,
        "verbose_mode": True,
        "sequencesToDomesticate": [
            {
                "sequence": "",
                "primerName": "",
                "mtkPartLeft": "",
                "mtkPartRight": ""
            }
        ]
    },
    "production": {
        "templateSequence": "",
        "species": "",
        "kozak": "",
        "max_mut_per_site": 1,
        "verbose_mode": False,
        "sequencesToDomesticate": []
    }
}
