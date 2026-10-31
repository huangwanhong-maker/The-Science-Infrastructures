"""OpenStreetMap's ``opening_hours`` read into the hours of a week.

The key's grammar is OpenStreetMap's specification
(https://wiki.openstreetmap.org/wiki/Key:opening_hours/specification, version of 2024). This module
reads a subset of it, the part that describes a typical week, and says how much of each value it
understood.

1. **Rules.** A value is a sequence of rules. A rule after ``;`` is a *normal* rule: for the days it
   names it replaces what earlier rules said. A rule after ``,`` (where the comma cannot continue a
   list of weekdays or of times) is an *additional* rule: it adds to the days it names. A rule after
   ``||`` is a *fallback* rule, which applies when no other matches; it is not a weekly hour and is
   left out (the value is then ``partial``).
2. **Selectors.** Weekdays ``Mo Tu We Th Fr Sa Su``, single, as ranges (``Mo-Fr``, and ``Fr-Mo``
   across the week's end) and lists (``Mo,We,Fr-Su``); ``PH``, the public holidays, kept apart from
   the week; times ``hh:mm-hh:mm`` and lists of them. A rule naming no weekday applies to every day; a
   rule naming weekdays and no time is open all day on them; ``24/7`` is open always.
3. **Past midnight.** A span whose end is not after its start runs into the next day (``22:00-02:00``
   ends at 26:00, 1560 minutes), as does an extended time (``18:00-26:00``); ``24:00`` is the day's
   end. A span of equal ends (``00:00-00:00``) is the whole day from its start. A day's hours are the
   spans that start on it: a later rule naming Tuesday replaces Tuesday's spans, not Monday's span
   that runs into Tuesday.
4. **States.** ``off`` and ``closed`` close the days named (in a normal or an additional rule);
   ``open`` is the default; ``unknown`` is left out. Comments (``"..."``) are ignored.
5. **What is left out.** Rules with a year, month, date, week or ``easter`` selector, with the n-th
   weekday of a month (``Su[1]``) or a day offset, with ``SH`` (school holidays) alone, and spans
   with a variable time (``sunrise``, ``sunset``, ``dawn``, ``dusk``) or an open end (``18:00+``) are
   left out, the other rules kept, and the value marked ``partial``; ``SH`` beside weekdays is
   dropped from the list; a span with an open end after it (``17:00-22:00+``) is kept to its stated
   end, and the value marked ``partial``. Periods within a span (``10:00-16:00/01:00``) keep the span and mark the
   value ``partial``. A rule that cannot be read is left out in the same way; a value of which no
   rule can be read is ``unparsed``.

The tolerance for common slips is stated: case is ignored in keywords and weekdays, whitespace is
free, a single-digit hour is read (``9:00``), en and em dashes and the wave dashes (``〜``, ``～``)
are read as ``-``, full-width characters are folded by NFKC, and a colon after a weekday selector is
skipped (``Mo-Fr: 10:00-18:00``). Nothing else is guessed: ``10am-6pm`` and ``Mon-Fri`` are not read.

The result is a dictionary: ``hours``, a list of ``[day, open, close]`` with day 0 Monday to 6
Sunday, and open and close in minutes from that day's 00:00 local time (close may exceed 1440 past
midnight), ordered and with overlapping spans of a day joined; ``ph``, the public holidays' spans
``[[open, close], ...]`` (``[]`` closed), or ``None`` when the value does not say; ``status``, one
of ``ok``, ``partial``, ``unparsed`` or ``none`` (no value); and ``dropped``, the rules or parts left
out, each as its text and the reason.

Usage::

    from opening_hours import parse
    parse("Mo-Fr 10:00-19:00; Sa 10:00-17:00; Su,PH off")
    # {'hours': [[0, 600, 1140], ..., [5, 600, 1020]], 'ph': [], 'status': 'ok', 'dropped': []}
"""

from __future__ import annotations

import re
import unicodedata

DAYS = ("mo", "tu", "we", "th", "fr", "sa", "su")
MONTHS = ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
VARIABLE = ("sunrise", "sunset", "dawn", "dusk")
STATES = ("open", "closed", "off", "unknown")
PH = 7                      # the public holidays, kept as an eighth day
DAY = 1440                  # minutes in a day

TOKEN = re.compile(r"""
      (?P<comment>"[^"]*")
    | (?P<always>24/7)
    | (?P<time>\d{1,2}:\d{2})
    | (?P<number>\d+)
    | (?P<word>[A-Za-z]+)
    | (?P<sym>[-,:\[\]+/()])
    """, re.VERBOSE)
WORDS = {**{d: "day" for d in DAYS}, **{m: "month" for m in MONTHS}, **{v: "variable" for v in VARIABLE},
         **{s: "state" for s in STATES}, "ph": "ph", "sh": "sh", "week": "wide", "easter": "wide",
         "day": "offset", "days": "offset"}


class Unreadable(ValueError):
    """A rule this subset of the grammar cannot read."""


# ---- text -----------------------------------------------------------------------------------------

def normalise(text: str) -> str:
    """Full-width characters folded (NFKC), dashes and wave dashes read as hyphens."""
    text = unicodedata.normalize("NFKC", text)
    for dash in ("‐", "‑", "‒", "–", "—", "−", "〜", "~"):
        text = text.replace(dash, "-")
    return text.strip()


def sequences(text: str) -> list[tuple[str, str]]:
    """The rule sequences of a value, each with the separator before it ('', ';' or '||'), split
    outside comments."""
    out, start, sep, quoted, i = [], 0, "", False, 0
    while i < len(text):
        c = text[i]
        if c == '"':
            quoted = not quoted
        elif not quoted and c == ";":
            out.append((sep, text[start:i]))
            start, sep = i + 1, ";"
        elif not quoted and text.startswith("||", i):
            out.append((sep, text[start:i]))
            start, sep = i + 2, "||"
            i += 1
        i += 1
    out.append((sep, text[start:]))
    return [(s, t.strip()) for s, t in out if t.strip()]


def tokens(text: str) -> list[tuple[str, str]]:
    """(kind, value) of each token of a rule sequence; words lowered and classified."""
    out, i = [], 0
    while i < len(text):
        if text[i].isspace():
            i += 1
            continue
        m = TOKEN.match(text, i)
        if not m:
            raise Unreadable(f"unreadable at {text[i:i + 12]!r}")
        kind, value = m.lastgroup, m.group(m.lastgroup)
        if kind == "word":
            value = value.lower()
            if value not in WORDS:
                raise Unreadable(f"unknown word {m.group(kind)!r}")
            kind = WORDS[value]
        out.append((kind, value))
        i = m.end()
    return out


def _minutes(hhmm: str) -> int:
    h, m = (int(x) for x in hhmm.split(":"))
    if m >= 60 or h > 48:
        raise Unreadable(f"no such time {hhmm}")
    return 60 * h + m


# ---- rules ----------------------------------------------------------------------------------------

class Rule:
    """One rule: its days (None for every day), its spans (None for the whole day), its state, and
    what was left out of it (drop: the whole rule; partial: parts, the rest kept)."""

    def __init__(self, additional: bool):
        self.additional = additional
        self.days: list[int] | None = None
        self.spans: list[tuple[int, int]] | None = None
        self.state = "open"
        self.drop: str | None = None
        self.partial: list[str] = []


class _Reader:
    def __init__(self, toks: list[tuple[str, str]]):
        self.t, self.i = toks, 0

    def peek(self, k: int = 0) -> tuple:
        j = self.i + k
        return self.t[j] if j < len(self.t) else (None, None)

    def take(self) -> tuple:
        tok = self.peek()
        self.i += 1
        return tok

    def done(self) -> bool:
        return self.i >= len(self.t)


def rules(toks: list[tuple[str, str]]) -> list[Rule]:
    """The rules of one sequence: the first, then any additional rules after commas."""
    r = _Reader(toks)
    out = [_rule(r, additional=False)]
    while not r.done():
        if r.peek() == ("sym", ",") and r.peek(1)[0] in ("day", "ph", "sh", "number", "month", "wide",
                                                         "time", "always", "state", "comment"):
            r.take()
            out.append(_rule(r, additional=True))
        else:
            raise Unreadable(f"unexpected {r.peek()[1]!r}")
    return out


def _rule(r: _Reader, additional: bool) -> Rule:
    rule = Rule(additional)
    # ---- wide-range selectors: years, months, dates, weeks, easter (the rule is left out) ----
    wide = []
    while r.peek()[0] in ("number", "month", "wide") or (
            wide and r.peek() in (("sym", "-"), ("sym", "/"), ("sym", ",")) and r.peek(1)[0] in ("number", "month", "wide")):
        wide.append(r.take()[1])
    if wide:
        rule.drop = f"a year, month, date or week selector ({' '.join(wide)})"
    # ---- weekday selector ----
    named, days = False, []
    while r.peek()[0] in ("day", "ph", "sh"):
        named = True
        kind, value = r.take()
        if kind == "ph":
            days.append(PH)
        elif kind == "sh":
            rule.partial.append("school holidays (SH)")
        else:
            first = last = DAYS.index(value)
            if r.peek() == ("sym", "-") and r.peek(1)[0] == "day":
                r.take()
                last = DAYS.index(r.take()[1])
            if r.peek() == ("sym", "["):                    # the n-th weekday of a month
                while not r.done() and r.peek() != ("sym", "]"):
                    r.take()
                if r.done():
                    raise Unreadable("an unclosed [")
                r.take()
                rule.drop = rule.drop or "the n-th weekday of a month"
            if r.peek()[0] == "sym" and r.peek()[1] in "+-" and r.peek(1)[0] == "number" and r.peek(2)[0] == "offset":
                r.take(), r.take(), r.take()
                rule.drop = rule.drop or "a day offset"
            days.extend((first + k) % 7 for k in range((last - first) % 7 + 1))
        if r.peek() == ("sym", ",") and r.peek(1)[0] in ("day", "ph", "sh"):
            r.take()
            continue
        break
    if named:
        if not days:
            rule.drop = rule.drop or "school holidays (SH) alone"
        rule.days = sorted(set(days))
        if r.peek() == ("sym", ":"):                        # a slip: "Mo-Fr: 10:00-18:00"
            r.take()
    # ---- time selector ----
    if r.peek()[0] == "always":
        r.take()
        rule.spans = [(0, DAY)]
    elif r.peek()[0] in ("time", "variable") or r.peek() == ("sym", "("):
        rule.spans = []
        while True:
            start, start_variable = _time_point(r)
            if r.peek() == ("sym", "+"):
                r.take()
                rule.partial.append("an open end (+)")
                end = None
            else:
                if r.peek() != ("sym", "-"):
                    raise Unreadable("a time without an end")
                r.take()
                end, end_variable = _time_point(r)
                if r.peek() == ("sym", "+"):            # 17:00-22:00+: open on after 22:00, for how long unsaid
                    r.take()
                    rule.partial.append("an open end after a span (kept to its stated end)")
            periodic = False
            if r.peek() == ("sym", "/"):
                r.take()
                if r.take()[0] not in ("time", "number"):
                    raise Unreadable("a period without a length")
                periodic = True
            if end is not None:
                if start_variable or end_variable:
                    rule.partial.append("a variable time (sunrise, sunset, dawn, dusk)")
                else:
                    if start >= DAY:
                        raise Unreadable("a span starting after 24:00")
                    rule.spans.append((start, end if end > start else end + DAY))
                    if periodic:
                        rule.partial.append("a period within a span (kept as the span)")
            if r.peek() == ("sym", ",") and (r.peek(1)[0] in ("time", "variable") or r.peek(1) == ("sym", "(")):
                r.take()
                continue
            break
        if not rule.spans:
            rule.drop = rule.drop or "no fixed time: " + ", ".join(rule.partial)
    # ---- state and comments ----
    if r.peek()[0] == "state":
        state = r.take()[1]
        rule.state = "closed" if state in ("closed", "off") else state
        if state == "unknown":
            rule.drop = rule.drop or "the state unknown"
    while r.peek()[0] == "comment":
        r.take()
    if not named and rule.spans is None and not wide and rule.state == "open":
        rule.drop = rule.drop or "a comment only"
    return rule


def _time_point(r: _Reader) -> tuple[int, bool]:
    """A time: hh:mm in minutes, or a variable time (sunrise, (sunset-01:00)), whose minutes are not
    known; returns (minutes, variable)."""
    kind, value = r.peek()
    if kind == "time":
        r.take()
        return _minutes(value), False
    if kind == "variable":
        r.take()
        return 0, True
    if (kind, value) == ("sym", "("):
        while not r.done() and r.peek() != ("sym", ")"):
            r.take()
        if r.done():
            raise Unreadable("an unclosed (")
        r.take()
        return 0, True
    raise Unreadable("a time expected")


# ---- the week ---------------------------------------------------------------------------------------

def parse(value: str | None) -> dict:
    """The weekly hours of an opening_hours value (see the module's docstring)."""
    if value is None or not str(value).strip():
        return {"hours": [], "ph": None, "status": "none", "dropped": []}
    text = normalise(str(value))
    week: list[list[tuple[int, int]]] = [[] for _ in range(8)]
    ph_stated, partial, read, dropped = False, False, 0, []
    for sep, piece in sequences(text):
        if sep == "||":
            dropped.append(f"{piece}: a fallback rule")
            partial = True
            continue
        try:
            found = rules(tokens(piece))
        except (Unreadable, ValueError) as error:
            dropped.append(f"{piece}: {error}")
            partial = True
            continue
        read += 1
        for rule in found:
            if rule.drop:
                dropped.append(f"{piece}: {rule.drop}")
                partial = True
                continue
            for note in rule.partial:
                dropped.append(f"{piece}: {note}")
                partial = True
            days = list(range(7)) if rule.days is None else rule.days
            spans = [(0, DAY)] if rule.spans is None else rule.spans
            for d in days:
                if rule.state == "closed":
                    week[d] = []
                elif rule.additional:
                    week[d] = week[d] + spans
                else:
                    week[d] = list(spans)
                ph_stated |= d == PH
    if not read:
        return {"hours": [], "ph": None, "status": "unparsed", "dropped": dropped}
    hours = [[d, a, b] for d in range(7) for a, b in merge(week[d])]
    ph = [[a, b] for a, b in merge(week[PH])] if ph_stated else None
    return {"hours": hours, "ph": ph, "status": "partial" if partial else "ok", "dropped": dropped}


def merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The spans of one day in order, overlapping or touching ones joined."""
    out: list[list[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def weekly_minutes(hours: list[list[int]]) -> int:
    """The minutes open in a week; spans past midnight counted whole, overlaps across days once."""
    covered = set()
    for d, a, b in hours:
        covered.update((d * DAY + m) % (7 * DAY) for m in range(a, b))
    return len(covered)


def is_open(hours: list[list[int]], day: int, minute: int) -> bool:
    """Whether the hours are open at a minute (0–1439) of a weekday (0 Monday), spans from the day
    before that run past midnight included."""
    for d, a, b in hours:
        if d == day and a <= minute < b:
            return True
        if (d + 1) % 7 == day and a <= minute + DAY < b:
            return True
    return False
