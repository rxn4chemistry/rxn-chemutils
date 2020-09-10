import urllib.parse

INTERNAL_RENDERING_SERVICE = 'http://depict.ava19.zc2.ibm.com:8080/depict/bow'
EXTERNAL_RENDERING_SERVICE = 'https://www.simolecule.com/cdkdepict/depict/bow'


def smiles_depict_url(smiles: str, format: str = 'svg', use_internal_service=True) -> str:
    """
    Generate the URL for the depiction of a SMILES string.

    Args:
        smiles: smiles string to depict
        format: 'svg', 'pdf', 'png', etc.
        use_internal_service: whether to use the service deployed on ZC2 (True)
            or the one available on the simolecule website (False).

    Returns:
        URL string
    """
    rendering_service = INTERNAL_RENDERING_SERVICE if use_internal_service else EXTERNAL_RENDERING_SERVICE
    params = {
        'smi': smiles,
        'zoom': '1.0',
        'abbr': 'on',
        'hdisp': 'bridgehead',
        'showtitle': 'false',
        'annotate': 'none'
    }
    params_str = urllib.parse.urlencode(params)
    return f'{rendering_service}/{format}?{params_str}'
