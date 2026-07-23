"""Streamlit web interface for dna_concat DNA assembly."""

import csv
import tempfile
from io import StringIO
from pathlib import Path

import streamlit as st
import yaml

from dna_concat.assembly import assemble
from dna_concat.io import read_csv as read_csv_file
from dna_concat.io import read_fasta
from dna_concat.models import AssemblySpec, NamedSequence, Slot, SlotBehavior

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(page_title="dna-concat — DNA Assembly", layout="wide")

# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------

_BEHAVIOR_OPTIONS = ["fixed", "zip", "product", "variable_stuffer"]


def _init_state():
    defaults = {
        "name_template": "",
        "allow_zip_trim": False,
        "slots": [],
        "next_slot_id": 0,
        "results": None,
        "result_name_template": None,
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


_init_state()


def _add_slot(
    name: str = "",
    behavior: str = "product",
    source_type: str = "Inline",
    inline_text: str = "",
    sequences: list[dict] | None = None,
    counterpart: str | None = None,
    truncate_side: str | None = None,
) -> None:
    """Append a new slot to session state with a stable ID."""
    slot_id = st.session_state["next_slot_id"]
    st.session_state["next_slot_id"] = slot_id + 1
    st.session_state["slots"].append(
        {
            "id": slot_id,
            "name": name,
            "behavior": behavior,
            "source_type": source_type,
            "inline_text": inline_text,
            "sequences": sequences or [],
            "counterpart": counterpart,
            "truncate_side": truncate_side,
        }
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_inline_text(text: str) -> list[NamedSequence]:
    """Parse inline text: one name<TAB or SPACE>sequence per line."""
    results = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise ValueError(
                f"Each line must be 'name<TAB>sequence'. Got: {line!r}"
            )
        results.append(NamedSequence(name=parts[0], sequence=parts[1]))
    return results


def _sequences_from_slot(slot_data: dict, prefix: str) -> list[NamedSequence]:
    """Extract NamedSequence list from a slot's session state data."""
    src = slot_data["source_type"]

    if src == "Inline":
        return _parse_inline_text(slot_data["inline_text"])

    if src == "Upload FASTA":
        uploaded = st.session_state.get(f"{prefix}_fasta_file")
        if uploaded is None:
            raise ValueError(f"Slot '{slot_data['name']}': no FASTA file uploaded")
        return slot_data.get("sequences", [])

    if src == "Upload CSV":
        uploaded = st.session_state.get(f"{prefix}_csv_file")
        if uploaded is None:
            raise ValueError(f"Slot '{slot_data['name']}': no CSV file uploaded")
        return slot_data.get("sequences", [])

    raise ValueError(f"Unknown source type: {src}")


def _build_results_table(
    constructs: list, name_template: str
) -> list[dict]:
    """Build a list of dicts for st.dataframe."""
    rows = []
    for c in constructs:
        row = {
            "name": c.format_name(name_template),
            "full_sequence": c.full_sequence,
            "length": len(c.full_sequence),
        }
        for slot_name in c.slot_order:
            ns = c.assignments[slot_name]
            row[f"{slot_name}"] = ns.name
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Sidebar — YAML config upload
# ---------------------------------------------------------------------------

st.sidebar.header("Load Configuration")
config_file = st.sidebar.file_uploader(
    "Upload a YAML config file",
    type=["yaml", "yml"],
    help="Same format as the CLI config. Populates slots and settings below.",
)

if config_file is not None:
    if st.sidebar.button("Apply Config"):
        raw = yaml.safe_load(config_file.getvalue())
        if not isinstance(raw, dict):
            st.sidebar.error("Config file must be a YAML mapping.")
        else:
            st.session_state["name_template"] = raw.get("name_template", "") or ""
            st.session_state["allow_zip_trim"] = bool(raw.get("allow_zip_trim", False))
            st.session_state["slots"] = []
            st.session_state["results"] = None

            file_ref_warnings = []
            for slot_def in raw.get("slots", []):
                name = slot_def.get("name", "")
                behavior = slot_def.get("behavior", "product")
                counterpart = slot_def.get("counterpart")
                truncate_side = slot_def.get("truncate_side")
                source = slot_def.get("source", {})

                if "inline" in source:
                    inline_entries = source["inline"]
                    inline_text = "\n".join(
                        f"{e['name']}\t{e['sequence']}" for e in inline_entries
                    )
                    _add_slot(
                        name=name,
                        behavior=behavior,
                        source_type="Inline",
                        inline_text=inline_text,
                        counterpart=counterpart,
                        truncate_side=truncate_side,
                    )
                elif "fasta" in source or "csv" in source:
                    file_ref_warnings.append(name)
                    _add_slot(
                        name=name,
                        behavior=behavior,
                        source_type="Inline",
                        counterpart=counterpart,
                        truncate_side=truncate_side,
                    )
                else:
                    _add_slot(
                        name=name,
                        behavior=behavior,
                        source_type="Inline",
                        counterpart=counterpart,
                        truncate_side=truncate_side,
                    )

            if file_ref_warnings:
                st.sidebar.warning(
                    f"Slots with file references were added empty "
                    f"(external file paths cannot be resolved in the browser): "
                    f"{', '.join(file_ref_warnings)}. "
                    f"Please re-add their sequences via file upload."
                )
            st.rerun()

st.sidebar.markdown("---")
st.sidebar.markdown(
    "Or configure settings manually below and build "
    "slots in the main panel."
)

# ---------------------------------------------------------------------------
# Title
# ---------------------------------------------------------------------------

st.title("dna-concat — DNA Assembly")
st.caption("Web interface for combinatorial DNA construct generation")

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

st.header("1. Settings")

col1, col2 = st.columns(2)

with col1:
    name_template_val = st.text_input(
        "Name template",
        value=st.session_state["name_template"],
        help=(
            "Python format string for construct names using slot names as keys. "
            'E.g. "{promoter}_{orf}". Leave blank for auto-generation.'
        ),
        key="name_template_input",
    )
with col2:
    allow_zip_trim_val = st.checkbox(
        "Allow zip trim",
        value=st.session_state["allow_zip_trim"],
        help="Allow zip slots to be trimmed when they exceed the product count.",
        key="allow_zip_trim_input",
    )

st.session_state["name_template"] = name_template_val
st.session_state["allow_zip_trim"] = allow_zip_trim_val

# ---------------------------------------------------------------------------
# Slot builder
# ---------------------------------------------------------------------------

st.header("2. Slots")

if st.button("Add Slot", key="add_slot"):
    _add_slot()
    st.rerun()

slots_to_keep = []
for idx, slot_data in enumerate(st.session_state["slots"]):
    sid = slot_data["id"]
    prefix = f"slot_{sid}"
    display_name = slot_data["name"] or f"(unnamed slot {idx + 1})"

    with st.expander(
        f"Slot {idx + 1}: {display_name} [{slot_data['behavior']}]",
        expanded=True,
    ):
        # Top row: name, behavior, remove
        c1, c2, c3 = st.columns([3, 2, 1])
        with c1:
            slot_data["name"] = st.text_input(
                "Slot name",
                value=slot_data["name"],
                key=f"{prefix}_name",
            )
        with c2:
            beh_idx = (
                _BEHAVIOR_OPTIONS.index(slot_data["behavior"])
                if slot_data["behavior"] in _BEHAVIOR_OPTIONS
                else 0
            )
            slot_data["behavior"] = st.selectbox(
                "Behavior",
                _BEHAVIOR_OPTIONS,
                index=beh_idx,
                key=f"{prefix}_behavior",
            )
        with c3:
            st.markdown("<br>", unsafe_allow_html=True)
            remove = st.button("Remove", key=f"{prefix}_remove")

        if remove:
            continue  # skip, don't add to slots_to_keep

        if slot_data["behavior"] == "fixed":
            st.info("Fixed slots require exactly 1 sequence.")
        elif slot_data["behavior"] == "variable_stuffer":
            st.info(
                "Variable stuffer slots use exactly 1 sequence and truncate it "
                "per construct so the stuffer plus its counterpart slot always "
                "equal the untruncated stuffer length."
            )
            other_names = [
                s["name"]
                for s in st.session_state["slots"]
                if s["id"] != sid and s["name"]
            ]
            sc1, sc2 = st.columns(2)
            with sc1:
                if other_names:
                    cp_idx = (
                        other_names.index(slot_data["counterpart"])
                        if slot_data.get("counterpart") in other_names
                        else 0
                    )
                    slot_data["counterpart"] = st.selectbox(
                        "Counterpart slot",
                        other_names,
                        index=cp_idx,
                        key=f"{prefix}_counterpart",
                        help="The slot whose length is subtracted from the stuffer.",
                    )
                else:
                    st.warning("Add another named slot to use as the counterpart.")
                    slot_data["counterpart"] = None
            with sc2:
                ts_options = ["left", "right"]
                ts_idx = (
                    ts_options.index(slot_data["truncate_side"])
                    if slot_data.get("truncate_side") in ts_options
                    else 1  # default: trim the right end
                )
                slot_data["truncate_side"] = st.selectbox(
                    "Truncate side",
                    ts_options,
                    index=ts_idx,
                    key=f"{prefix}_truncate_side",
                    help="Which end of the stuffer sequence to trim.",
                )

        # Source type
        source_options = ["Inline", "Upload FASTA", "Upload CSV"]
        src_idx = (
            source_options.index(slot_data["source_type"])
            if slot_data["source_type"] in source_options
            else 0
        )
        slot_data["source_type"] = st.radio(
            "Source",
            source_options,
            index=src_idx,
            key=f"{prefix}_source",
            horizontal=True,
        )

        if slot_data["source_type"] == "Inline":
            slot_data["inline_text"] = st.text_area(
                "Sequences (one per line: name<TAB>sequence)",
                value=slot_data.get("inline_text", ""),
                height=120,
                key=f"{prefix}_inline",
                help="Enter one sequence per line: name followed by tab or space, then the DNA sequence.",
            )

        elif slot_data["source_type"] == "Upload FASTA":
            fasta_file = st.file_uploader(
                "Upload FASTA",
                type=["fa", "fasta", "fna"],
                key=f"{prefix}_fasta_file",
            )
            if fasta_file is not None:
                suffix = Path(fasta_file.name).suffix
                tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
                tmp.write(fasta_file.getvalue())
                tmp.close()
                try:
                    seqs = read_fasta(tmp.name)
                    slot_data["sequences"] = seqs
                    st.success(f"{len(seqs)} sequence(s) loaded from FASTA.")
                except Exception as e:
                    st.error(f"Error reading FASTA: {e}")
                    slot_data["sequences"] = []
                finally:
                    Path(tmp.name).unlink(missing_ok=True)

        elif slot_data["source_type"] == "Upload CSV":
            csv_file = st.file_uploader(
                "Upload CSV",
                type=["csv", "tsv", "txt"],
                key=f"{prefix}_csv_file",
            )
            cc1, cc2, cc3 = st.columns(3)
            with cc1:
                name_col = st.text_input(
                    "Name column",
                    value="name",
                    key=f"{prefix}_csv_name_col",
                )
            with cc2:
                seq_col = st.text_input(
                    "Sequence column",
                    value="sequence",
                    key=f"{prefix}_csv_seq_col",
                )
            with cc3:
                delimiter = st.text_input(
                    "Delimiter",
                    value=",",
                    key=f"{prefix}_csv_delim",
                )
            if csv_file is not None:
                suffix = Path(csv_file.name).suffix
                tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
                tmp.write(csv_file.getvalue())
                tmp.close()
                try:
                    seqs = read_csv_file(
                        tmp.name,
                        name_column=name_col,
                        sequence_column=seq_col,
                        delimiter=delimiter,
                    )
                    slot_data["sequences"] = seqs
                    st.success(f"{len(seqs)} sequence(s) loaded from CSV.")
                except Exception as e:
                    st.error(f"Error reading CSV: {e}")
                    slot_data["sequences"] = []
                finally:
                    Path(tmp.name).unlink(missing_ok=True)

        slots_to_keep.append(slot_data)

# Update session state if a slot was removed
if len(slots_to_keep) != len(st.session_state["slots"]):
    st.session_state["slots"] = slots_to_keep
    st.session_state["results"] = None
    st.rerun()
else:
    st.session_state["slots"] = slots_to_keep

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

st.header("3. Assemble")

dry_run = st.checkbox(
    "Dry run",
    help="Show slot summary and expected construct count without generating full results.",
    key="dry_run_check",
)

if st.button("Run Assembly", type="primary", key="run_assembly"):
    if not st.session_state["slots"]:
        st.error("Please add at least one slot.")
        st.stop()

    # Build Slot objects
    built_slots = []
    for slot_data in st.session_state["slots"]:
        if not slot_data["name"]:
            st.error("All slots must have a name.")
            st.stop()

        try:
            if slot_data["source_type"] == "Inline":
                seqs = _parse_inline_text(slot_data.get("inline_text", ""))
            else:
                seqs = slot_data.get("sequences", [])

            if not seqs:
                st.error(
                    f"Slot '{slot_data['name']}' has no sequences. "
                    f"Please add sequences or upload a file."
                )
                st.stop()

            behavior = SlotBehavior(slot_data["behavior"])
            extra = {}
            if behavior is SlotBehavior.VARIABLE_STUFFER:
                extra["counterpart"] = slot_data.get("counterpart") or None
                extra["truncate_side"] = slot_data.get("truncate_side") or None
            built_slots.append(
                Slot(name=slot_data["name"], behavior=behavior, sequences=seqs, **extra)
            )
        except ValueError as e:
            st.error(f"Slot '{slot_data['name']}': {e}")
            st.stop()

    # Build AssemblySpec
    template = st.session_state["name_template"] or None
    try:
        spec = AssemblySpec(
            slots=built_slots,
            name_template=template,
            allow_zip_trim=st.session_state["allow_zip_trim"],
        )
    except ValueError as e:
        st.error(f"Assembly configuration error: {e}")
        st.stop()

    # Run assembly
    try:
        constructs = assemble(spec)
    except ValueError as e:
        st.error(f"Assembly error: {e}")
        st.stop()

    if dry_run:
        st.subheader("Dry Run Summary")
        st.markdown(f"**Slots:** {len(spec.slots)}")
        for s in spec.slots:
            st.markdown(
                f"- **{s.name}** ({s.behavior.value}): "
                f"{len(s.sequences)} sequence(s)"
            )
        st.markdown(f"**Total constructs:** {len(constructs)}")
        st.markdown(f"**Name template:** `{spec.name_template}`")
        if constructs:
            first = constructs[0]
            st.markdown(
                f"**First construct:** {first.format_name(spec.name_template)} "
                f"({len(first.full_sequence)} bp)"
            )
    else:
        st.session_state["results"] = constructs
        st.session_state["result_name_template"] = spec.name_template
        st.rerun()

# ---------------------------------------------------------------------------
# Results display
# ---------------------------------------------------------------------------

if st.session_state["results"] is not None:
    constructs = st.session_state["results"]
    name_tmpl = st.session_state["result_name_template"] or "construct"

    st.header("Results")
    st.success(f"{len(constructs)} construct(s) generated.")

    # Results table
    if constructs:
        rows = _build_results_table(constructs, name_tmpl)
        st.dataframe(rows, use_container_width=True)

    # Downloads
    st.subheader("Download Results")
    dl_col1, dl_col2 = st.columns(2)

    # FASTA download
    fasta_buf = StringIO()
    for c in constructs:
        name = c.format_name(name_tmpl)
        fasta_buf.write(f">{name}\n{c.full_sequence}\n")
    with dl_col1:
        st.download_button(
            "Download FASTA",
            data=fasta_buf.getvalue(),
            file_name="constructs.fasta",
            mime="text/plain",
        )

    # CSV download
    if constructs:
        csv_buf = StringIO()
        slot_names = constructs[0].slot_order
        fieldnames = ["name", "full_sequence"]
        for sn in slot_names:
            fieldnames.extend([f"{sn}_name", f"{sn}_seq"])

        writer = csv.DictWriter(csv_buf, fieldnames=fieldnames)
        writer.writeheader()
        for c in constructs:
            row = {
                "name": c.format_name(name_tmpl),
                "full_sequence": c.full_sequence,
            }
            for sn in slot_names:
                ns = c.assignments[sn]
                row[f"{sn}_name"] = ns.name
                row[f"{sn}_seq"] = ns.sequence
            writer.writerow(row)

        with dl_col2:
            st.download_button(
                "Download CSV",
                data=csv_buf.getvalue(),
                file_name="constructs.csv",
                mime="text/csv",
            )
