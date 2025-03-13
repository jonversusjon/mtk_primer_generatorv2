# config/test_config.py

CONFIG = {
    "development": {
        "templateSequence": "",
        "species": "",
        "kozak": "",
        "max_mut_per_site": 1,
        "verbose_mode": True,
        "sequencesToDomesticate": [
            {
                "sequence": "",
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
