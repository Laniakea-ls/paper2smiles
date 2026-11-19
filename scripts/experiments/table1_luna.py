"""Collect independent vision candidates for the Table 1 experiment."""


def main():
    """Run the local source-review experiment."""
    import getpass
    from pathlib import Path

    from paper2smiles.extraction.core import save
    from paper2smiles.extraction.model import Luna

    root = Path("runs/gpr88_table1_reconstruction").resolve()
    llm = Luna(
        getpass.getpass("API key (hidden): "),
        "https://newapi.livzon.cn/v1",
        "deepseek-v4.1-flash",
        root,
    )
    prompt = (
        "Read ONLY the supplied original Table 1 and source PDF page. Test "
        "extraction of the generic scaffold plus each R1 fragment for "
        "compounds 4–9. Do not consult other literature or invent names. "
        'Return JSON {"scaffold":{"smiles":"SMILES with exactly one [*:1] at '
        'the aromatic R1 bond","stereo_evidence":"what is actually drawn and '
        'corroborated in text","uncertainties":[]},"rows":[{"compound_id":"4",'
        '"r1_smiles":"fragment with exactly one [*:1] representing the bond '
        'back to the aromatic scaffold","connectivity_description":"explicit '
        'atoms and attachment orientation","stereo_evidence":"which '
        'stereocenters are specified or unspecified/racemic","ec50_printed":".'
        '..","uncertainties":[]}],"warnings":[]}. The wavy slash crossing the '
        "attachment bond marks attachment, do not treat it as an extra atom. "
        "Preserve all explicit stereochemistry. Do not invent absolute "
        "stereochemistry for unspecified centers; record mixture information "
        "separately. The generic scaffold includes the entire fixed amide "
        "cap. Return a tentative candidate, not a claim of verified "
        "correctness."
    )
    r = llm.ask(
        "table1_independent_vision",
        prompt,
        [root / "table1.png", root / "scaffold.png", root / "page2.png"],
        max_tokens=6000,
    )
    save(root / "luna_candidates.json", r)
    print("Saved independent Luna candidate")


if __name__ == "__main__":
    main()
