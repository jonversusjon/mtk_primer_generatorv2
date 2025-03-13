from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent

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
                "primerName": "",
                "mtkPartLeft": "",
                "mtkPartRight": ""
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
