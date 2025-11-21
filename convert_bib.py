import argparse
import re
from pathlib import Path

MONTHS = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]

JOURNAL_ABBREVIATIONS = {
    "IEEE Transactions on Wireless Communications": "IEEE_J_WCOM",
    "IEEE Transactions on Communications": "IEEE_J_COM",
    "IEEE Transactions on Signal Processing": "IEEE_J_SP",
    "IEEE Transactions on Information Theory": "IEEE_J_IT",
}


STOPWORDS = {
    "a",
    "an",
    "and",
    "as",
    "at",
    "for",
    "from",
    "in",
    "of",
    "on",
    "or",
    "the",
    "to",
}


def wrap_acronyms(title: str) -> str:
    def transform(token: str, position: int) -> str:
        # Remove surrounding punctuation for acronym detection
        stripped = token.strip("{}(),.;:")
        if stripped.isupper() and len(stripped) > 1:
            # Preserve original token punctuation while wrapping acronym
            return token.replace(stripped, f"{{{stripped}}}")

        # Title-case normal words while respecting existing hyphens
        parts = token.split("-")
        titled = "-".join(p[:1].upper() + p[1:].lower() if p else "" for p in parts)
        lower_candidate = titled.lower()
        if position > 0 and lower_candidate in STOPWORDS:
            return lower_candidate
        return titled

    return " ".join(transform(tok, idx) for idx, tok in enumerate(title.split()))


def detect_month(number: str | None) -> str | None:
    if number is None:
        return None
    try:
        idx = int(number)
    except ValueError:
        return None
    if 1 <= idx <= 12:
        return MONTHS[idx - 1]
    return None


def generate_key(author_field: str, year: str, title: str) -> str:
    first_author = author_field.split(" and ")[0]
    surname = first_author.split(",")[0].strip()
    first_word = title.split()[0].lower() if title else ""
    return f"{surname}{year}{first_word}"


def parse_entry(text: str) -> tuple[str, dict[str, str]]:
    header_match = re.match(r"@(\w+)\s*\{\s*([^,]+)", text)
    if not header_match:
        raise ValueError("Unsupported entry format")
    entry_type, key = header_match.groups()
    fields: dict[str, str] = {}
    for field, value in re.findall(r"(\w+)\s*=\s*\{([^{}]*)\}", text):
        fields[field.lower()] = value.strip()
    return entry_type, key, fields


def format_entry(entry_type: str, key: str, fields: dict[str, str]) -> str:
    order = [
        "author",
        "journal",
        "title",
        "year",
        "month",
        "volume",
        "number",
        "pages",
        "keywords",
        "doi",
    ]
    lines = [f"@{entry_type.upper()}{{{key},"]
    for field in order:
        if field in fields:
            lines.append(f"  {field}={{{fields[field]}}},")
    # remove trailing comma on last line
    if lines:
        lines[-1] = lines[-1].rstrip(',')
    lines.append("}")
    return "\n".join(lines)


def convert_entry(text: str) -> str:
    entry_type, key, fields = parse_entry(text)

    # Standardize journal
    journal = fields.get("journal")
    if journal in JOURNAL_ABBREVIATIONS:
        fields["journal"] = JOURNAL_ABBREVIATIONS[journal]

    # Month from number if applicable
    month = detect_month(fields.get("number"))
    if month:
        fields["month"] = month

    # Title casing and acronym preservation
    if "title" in fields:
        fields["title"] = wrap_acronyms(fields["title"])

    # Generate new key
    if "author" in fields and "year" in fields and "title" in fields:
        key = generate_key(fields["author"], fields["year"], fields["title"])

    return format_entry(entry_type, key, fields)


def split_entries(content: str) -> list[str]:
    # Simple splitter assuming entries are separated by blank lines or newlines
    matches = re.findall(r"@\w+\s*\{[^@]*\}", content, re.DOTALL)
    return matches


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert IEEE bib entries to abbreviated format.")
    parser.add_argument("input", type=Path, help="Path to input .bib file")
    parser.add_argument("--output", type=Path, help="Optional output file; defaults to stdout")
    args = parser.parse_args()

    content = args.input.read_text(encoding="utf-8")
    entries = split_entries(content)
    if not entries:
        raise SystemExit("No bib entries found in input file")

    converted = "\n\n".join(convert_entry(entry) for entry in entries)

    if args.output:
        args.output.write_text(converted + "\n", encoding="utf-8")
    else:
        print(converted)


if __name__ == "__main__":
    main()
