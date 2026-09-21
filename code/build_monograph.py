# -*- coding: utf-8 -*-
r"""
build_monograph.py -- assemble the ANCCFT papers into one LaTeX book.

    python build_monograph.py analyse          # report macro / bibkey conflicts
    python build_monograph.py build <outdir>   # write <outdir>/main.tex, chapters/, macros, bib

Each paper (amsart or article) is parsed into title, abstract, preamble macros, body
(between \begin{document} and \begin{thebibliography}, minus \maketitle / abstract) and
\bibitem entries.  The book has one preamble (union of macros, first definition wins;
chapters whose own definition differs re-declare it locally), one family of theorem
environments numbered chapter.section.k (so "Thm 3.7 of Paper I" is Theorem 1.3.7 of
Chapter 1), and one merged bibliography.  Inter-paper references written as literal
"[V, Thm.~3.1]" are rewritten to "[Ch.~\ref{chap:V}, Thm.~3.1]".
"""
import os, re, sys, collections, shutil

ROOT = os.path.dirname(os.path.abspath(__file__))

# The same script runs in two layouts: beside the working folders ("Paper 1", "Phys 1", ...)
# and inside the repository (code/build_monograph.py with papers/01, applications/...).
REPO_MAP = {f"Paper {i}": f"papers/{i:02d}" for i in range(1, 21)}
REPO_MAP.update({"Phys 1": "applications/phys-bloch-mass", "Holo 1": "applications/holo-boundary-blind",
                 "Net 1": "applications/net-tower-lambda", "Ctrl 1": "applications/ctrl-hidden-eigenvalues"})


def locate(folder):
    """the directory holding a paper's sources, in either layout"""
    for cand in (os.path.join(ROOT, folder),
                 os.path.join(ROOT, REPO_MAP.get(folder, folder)),
                 os.path.join(os.path.dirname(ROOT), REPO_MAP.get(folder, folder))):
        if os.path.isdir(cand):
            return cand
    raise FileNotFoundError(f"no source folder for {folder!r} relative to {ROOT}")


# (part title, [(folder, tex, roman/label)])  -- the monograph plan
# Volume I: the twenty papers and the four non-biological applied notes.  The Bio strand
# (Bio 1--11) is Volume II, "Algebraic Genomics", built separately (see build_volume2 below).
PLAN = [
    ("The bridge over the tree: two channels and one obstruction", [
        ("Paper 1", "anccft.tex", "I"), ("Paper 2", "anccft2.tex", "II"),
        ("Paper 6", "anccft6.tex", "VI"), ("Paper 3", "anccft3.tex", "III")]),
    ("Towers: commutative Iwasawa theory and its failure", [
        ("Paper 4", "anccft4.tex", "IV"), ("Paper 9", "anccft9.tex", "IX"),
        ("Paper 14", "anccft14.tex", "XIV"), ("Paper 19", "anccft19.tex", "XIX"),
        ("Paper 5", "anccft5.tex", "V"), ("Paper 7", "anccft7.tex", "VII"),
        ("Paper 13", "anccft13.tex", "XIII")]),
    ("Operator algebras of the tower", [
        ("Paper 8", "anccft8.tex", "VIII"), ("Paper 17", "anccft17.tex", "XVII"),
        ("Paper 16", "anccft16.tex", "XVI")]),
    ("Higher rank: $\\widetilde A_2$ and $\\widetilde A_d$ buildings", [
        ("Paper 10", "anccft10.tex", "X"), ("Paper 11", "anccft11.tex", "XI"),
        ("Paper 12", "anccft12.tex", "XII"), ("Paper 15", "anccft15.tex", "XV"),
        ("Paper 18", "anccft18.tex", "XVIII"), ("Paper 20", "anccft20.tex", "XX")]),
    ("Applications to physics, networks and control", [
        ("Phys 1", "bloch_mass.tex", "Phys1"), ("Holo 1", "boundary_blind.tex", "Holo1"),
        ("Net 1", "tower_lambda.tex", "Net1"), ("Ctrl 1", "hidden_eigs.tex", "Ctrl1")]),
]
ORDER = [(f, t, l) for _, chs in PLAN for f, t, l in chs]

MACRO_RE = re.compile(r"\\(?:re)?newcommand\*?\{?\\(\w+)\}?(\[\d\])?(\[[^\]]*\])?\{")
MATHOP_RE = re.compile(r"\\DeclareMathOperator\*?\{\\(\w+)\}\{([^}]*)\}")
THM_RE = re.compile(r"\\newtheorem\*?\{(\w+)\}(?:\[(\w+)\])?\{([^}]*)\}(?:\[(\w+)\])?")
ROMAN = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII",
         "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"}


def balanced(s, i):
    depth, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return j + 1
        j += 1
    raise ValueError("unbalanced")


def parse(path):
    src = open(path, encoding="utf-8").read()
    pre, _, rest = src.partition("\\begin{document}")
    body, _, tail = rest.partition("\\begin{thebibliography}")
    m = re.search(r"\\title(\[[^\]]*\])?\s*\{", src)
    title = src[m.end():balanced(src, m.end() - 1) - 1] if m else os.path.basename(path).replace("_", "\\_")
    for cmd in ("title", "author", "date", "thanks", "address", "email", "subjclass", "keywords", "curraddr",
                "urladdr", "affiliation", "altaffiliation", "preprint", "pacs", "orcid", "homepage"):
        while True:
            mm = re.search(r"\\" + cmd + r"(\[[^\]]*\])?\s*\{", body)
            if not mm:
                break
            body = body[:mm.start()] + body[balanced(body, mm.end() - 1):]
    title = re.sub(r"\\\\\s*", " ", title).replace("\\bfseries", "").strip()
    macros = collections.OrderedDict()
    for mm in MACRO_RE.finditer(pre):
        macros[mm.group(1)] = pre[mm.start():balanced(pre, mm.end() - 1)]
    for mm in MATHOP_RE.finditer(pre):
        macros[mm.group(1)] = f"\\newcommand{{\\{mm.group(1)}}}{{\\operatorname{{{mm.group(2)}}}}}"
    thms = [(t.group(1), t.group(3)) for t in THM_RE.finditer(pre)]
    am = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", body, re.S)
    abstract = am.group(1).strip() if am else ""
    body = re.sub(r"\\begin\{abstract\}.*?\\end\{abstract\}", "", body, flags=re.S)
    body = re.sub(r"\\maketitle|\\tableofcontents", "", body)
    body = body.replace("\\end{document}", "")
    # a paper that sets its bibliography in \begingroup\small ... \endgroup leaves the opener
    # behind once the bibliography is cut off; drop it (the closer went with the bibliography)
    body = re.sub(r"\\begingroup\s*(\\(?:small|footnotesize|scriptsize))?\s*$", "", body.rstrip())
    assert body.count("\\begingroup") == body.count("\\endgroup"), path
    bib = collections.OrderedDict()
    if tail:
        tail = tail.split("\\end{thebibliography}")[0]
        parts = re.split(r"\\bibitem\{([^}]*)\}", tail)
        for k in range(1, len(parts), 2):
            bib[parts[k]] = " ".join(parts[k + 1].split())
    cites = set(re.findall(r"\\cite[tp]?\*?(?:\[[^\]]*\])?\{([^}]*)\}", body))
    return dict(title=title, macros=macros, thms=thms, abstract=abstract, body=body, bib=bib, cites=cites)


def norm(text):
    return " ".join(text.split())


def analyse(papers):
    print("=== titles ===")
    for lab, p in papers.items():
        print(f"  {lab:6s} {p['title'][:90]}   ({len(p['body'])//1000} kB body, {len(p['bib'])} refs)")
    print("\n=== macro conflicts ===")
    defs = collections.defaultdict(dict)
    for lab, p in papers.items():
        for name, text in p["macros"].items():
            defs[name].setdefault(norm(text), []).append(lab)
    for name, variants in sorted(defs.items()):
        if len(variants) > 1:
            print(f"  \\{name}: " + " | ".join(f"{t[:50]} <- {','.join(l)}" for t, l in variants.items()))
    print("\n=== theorem environments ===", dict(collections.Counter(e for p in papers.values() for e, _ in p["thms"])))
    keys = collections.defaultdict(dict)
    for lab, p in papers.items():
        for k, t in p["bib"].items():
            keys[k].setdefault(t, []).append(lab)
    print("\n=== bibliography key clashes ===")
    for k, variants in sorted(keys.items()):
        if len(variants) > 1:
            print(f"  {k}: " + " | ".join(f"{t[:40]} <- {','.join(l)}" for t, l in variants.items()))


# self-citations of unpublished parts of this corpus become chapter references.  The papers
# key them three ways: by bare roman numeral (\cite{II}), by ANCFT<n> and by Paper<n>; the
# applied notes by their own labels.  Every such key is dropped from the merged bibliography.
ROMANS = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII",
          "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"]
SELF = {r: r for r in ROMANS}
SELF.update({f"ANCFT{i + 1}": r for i, r in enumerate(ROMANS)})
SELF.update({f"Paper{i + 1}": r for i, r in enumerate(ROMANS)})
SELF.update({"Phys1": "Phys1", "Net1": "Net1", "Holo1": "Holo1", "Ctrl1": "Ctrl1"})
SELF.update({f"Bio{i}": f"Bio{i}" for i in range(1, 12)})
SELF_SERIES = "ANCFT"      # cited as \cite[V]{ANCFT}, \cite[X--XII]{ANCFT}


def chapter_refs(labels):
    refs = [f"\\ref{{chap:{l}}}" for l in labels]
    return ("Ch.~" if len(refs) == 1 else "Chs.~") + ", ".join(refs)


def rewrite_cites(body, mapping):
    def rep(m):
        cmd, opt, keys = m.group(1), m.group(2), [k.strip() for k in m.group(3).split(",")]
        keys = [mapping.get(k, k) for k in keys]
        ext = [k for k in keys if k not in SELF and k != SELF_SERIES]
        selfk = [k for k in keys if k in SELF or k == SELF_SERIES]
        if not selfk:
            return m.group(0)
        opt_text = opt[1:-1] if opt else ""
        if SELF_SERIES in selfk and opt_text:
            # optional argument names the papers: "V", "X--XII", "I, III"
            parts = []
            for tok in re.split(r",\s*", opt_text):
                rng = re.fullmatch(r"([IVX]+)\s*-{1,3}\s*([IVX]+)", tok.strip())
                if rng and rng.group(1) in ROMAN and rng.group(2) in ROMAN:
                    parts.append(f"\\ref{{chap:{rng.group(1)}}}--\\ref{{chap:{rng.group(2)}}}")
                elif tok.strip() in ROMAN:
                    parts.append(f"\\ref{{chap:{tok.strip()}}}")
                else:
                    parts.append(tok.strip())
            out = "[" + ("Ch.~" if len(parts) == 1 and "--" not in parts[0] else "Chs.~") + ", ".join(parts) + "]"
        else:
            labels = [SELF[k] for k in selfk if k != SELF_SERIES]
            if SELF_SERIES in selfk:
                labels = ["I", "XVI"]   # the series as a whole: Chapters I--XVI
                inner = f"Chs.~\\ref{{chap:I}}--\\ref{{chap:XVI}}"
            else:
                inner = chapter_refs(labels)
            out = "[" + inner + (", " + opt_text if opt_text else "") + "]"
        if ext:
            out = cmd + "{" + ",".join(ext) + "} and " + out
        return out
    return re.sub(r"(\\(?:online)?cite[tp]?\*?)(\[[^\]]*\])?\{([^}]*)\}", rep, body)


def _chref(r):
    return "\\ref{chap:" + r + "}"


# --- PDF bookmarks: hyperref cannot put math into a bookmark string, so every chapter, part
# and section title that contains math or a macro is wrapped as \texorpdfstring{tex}{plain}.
PDF_REP = {
    r"\widetilde A_2": "A2", r"\Atwo": "A2", r"\widetilde A_d": "Ad", r"\Atilde d": "Ad",
    r"\Atilde{d}": "Ad", r"\Atilde{2}": "A2", r"\Atilde 2": "A2", r"\Atilde{3}": "A3",
    r"\mathcal{L}": "L", r"\mathcal L": "L", r"\cL": "L", r"\mathrm{PGL}": "PGL", r"\PGL": "PGL",
    r"\mathrm{GL}": "GL", r"\GL": "GL", r"\mathrm{Sp}": "Sp", r"\PPone": "P1", r"\mathbb{P}": "P",
    r"\Zl": "Zl", r"\Zp": "Zp", r"\Qp": "Qp", r"\Ql": "Ql", r"\KK": "KK", r"\Z": "Z", r"\Q": "Q",
    r"\K": "K", r"\lambda": "lambda", r"\eta": "eta", r"\ell": "l", r"\Delta": "Delta",
    r"\Theta": "Theta", r"\Psi": "Psi", r"\tau": "tau", r"\alpha": "alpha", r"\sigma": "sigma",
    r"\kappa": "kappa", r"\mu": "mu", r"\nu": "nu", r"\infty": "inf", r"\rtimes": " x ",
    r"\otimes": " x ", r"\oplus": " + ", r"\to": " -> ", r"\emph": "", r"\textbf": "", r"\textit": "",
    r"\texttt": "", r"\mathrm": "", r"\operatorname": "", r"\ ": " ", r"\,": " ", r"\;": " ",
    r"\&": "&", r"\%": "%", r"\_": "_", "~": " ", "--": "-",
}


def pdf_plain(s):
    """an ASCII rendering of a title for the PDF bookmark"""
    s = re.sub(r"\\texorpdfstring\{((?:[^{}]|\{[^{}]*\})*)\}\{((?:[^{}]|\{[^{}]*\})*)\}", r"\2", s)
    for k in sorted(PDF_REP, key=len, reverse=True):
        s = s.replace(k, PDF_REP[k])
    s = re.sub(r"\\[A-Za-z]+\*?", "", s)
    s = re.sub(r"[\$\{\}\^]", "", s).replace("_", "")
    return re.sub(r"\s+", " ", s).strip()


def pdfsafe(title):
    if "\\texorpdfstring" in title or ("$" not in title and "\\" not in title):
        return title
    return "\\texorpdfstring{" + title + "}{" + pdf_plain(title) + "}"


def wrap_section_titles(body):
    r"""apply pdfsafe to the argument of every \section / \section* in a chapter body"""
    out, i = [], 0
    for m in re.finditer(r"\\section\*?\{", body):
        if m.start() < i:
            continue
        end = balanced(body, m.end() - 1)
        arg = body[m.end():end - 1]
        out.append(body[i:m.end()])
        out.append(pdfsafe(arg))
        out.append("}")
        i = end
    out.append(body[i:])
    return "".join(out)


def rewrite_paper_refs(body):
    r"""Turn the papers' own cross-references, written as roman numerals, into chapter
    references.  Only unambiguous citation shapes are rewritten: a BARE [X] can be a class
    in a module (\partial[\Theta]=[X] in the chapter from Paper XIII, and [x], [y] in the
    one from Paper VIII), so it is touched only when followed by punctuation, a space or an
    apostrophe, and never when preceded by a word character, a $ or a brace."""
    R = "|".join(sorted(ROMAN, key=len, reverse=True))   # longest first: XIII before XI before I
    # ranges, including the half-rewritten ones an earlier pass could leave
    body = re.sub(r"\[(" + R + r")\]--\[(" + R + r")\]",
                  lambda m: "Chs.~" + _chref(m.group(1)) + "--" + _chref(m.group(2)), body)
    body = re.sub(r"\[(" + R + r")\]--\[Ch\.~\\ref\{chap:(" + R + r")\}\]",
                  lambda m: "Chs.~" + _chref(m.group(1)) + "--" + _chref(m.group(2)), body)
    # two-paper lists: [I, III]
    body = re.sub(r"\[(" + R + r"),\s*(" + R + r")\]",
                  lambda m: "[Chs.~" + _chref(m.group(1)) + ", " + _chref(m.group(2)) + "]", body)
    # anything after a comma is a locator: [V, Thm.~3.1], [II,~1.4], [III, proof of Thm.~2.1]
    body = re.sub(r"\[(" + R + r"),\s*", lambda m: "[Ch.~" + _chref(m.group(1)) + ", ", body)
    # {\rm[I]} -- the \rm wrapper marks a citation
    body = re.sub(r"\\rm\s*\[(" + R + r")\]",
                  lambda m: "\\rm[Ch.~" + _chref(m.group(1)) + "]", body)
    # bare [V], and the possessive [X]'s, only in safe surroundings
    body = re.sub(r"(?<![\w$^_{\\])\[(" + R + r")\](?=['\s.,;:)])",
                  lambda m: "[Ch.~" + _chref(m.group(1)) + "]", body)
    return body


MAIN_TEMPLATE = r"""% main.tex -- generated by build_monograph.py; edit frontmatter/*.tex, not chapters/*.tex
\documentclass[11pt,twoside,openright]{book}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage[margin=1.05in]{geometry}
\usepackage{booktabs,array,longtable}
\usepackage{graphicx}
\usepackage{microtype}
\usepackage{url}
\usepackage[hidelinks]{hyperref}
\raggedbottom
\emergencystretch=2em
\setcounter{tocdepth}{1}
\setcounter{secnumdepth}{2}

\numberwithin{equation}{section}
\theoremstyle{plain}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{conjecture}[theorem]{Conjecture}
\newtheorem{hypothesis}[theorem]{Hypothesis}
\newtheorem{problem}[theorem]{Problem}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{algorithm}[theorem]{Algorithm}
\newtheorem{construction}[theorem]{Construction}
\newtheorem{convention}[theorem]{Convention}
\newtheorem{example}[theorem]{Example}
\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}

\newenvironment{chapterabstract}{\begin{quotation}\small\noindent\textbf{Summary.}\ }{\end{quotation}\medskip}
% compatibility layer for the revtex-style applied note
\usepackage{bm}
\providecommand{\onlinecite}[1]{\cite{#1}}
\newenvironment{ruledtabular}{}{}
\providecommand{\colrule}{\midrule}
\newenvironment{acknowledgments}{\section*{Acknowledgments}}{}
\input{macros}

\title{\bfseries Algorithmic Non-commutative Class Field Theory\\[6pt]
\large Hecke operators, sandpile groups and operator algebras on trees and buildings,\\ with applications}
\author{Ruqing Chen\\[4pt]\normalsize GUT Geoservice Inc., Rivi\`ere-Beaudette, Qu\'ebec, Canada\\
\normalsize\texttt{ruqing@hotmail.com}}
\date{September 2026\\[10pt]\normalsize Sources, code, data and figures:\\
\url{https://github.com/Ruqing1963/anccft}}

\begin{document}
\frontmatter
\maketitle
\input{frontmatter/preface}
\tableofcontents
\input{frontmatter/concordance}
\mainmatter
\input{frontmatter/introduction}
%%PARTS%%
\input{frontmatter/state}
\backmatter
\input{bibliography}
\end{document}
"""


def build(papers, outdir):
    os.makedirs(os.path.join(outdir, "chapters"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "frontmatter"), exist_ok=True)
    # unified macros
    macros = collections.OrderedDict()
    for lab, p in papers.items():
        for name, text in p["macros"].items():
            macros.setdefault(name, text)
    # merged bibliography
    bib, remap = collections.OrderedDict(), {}
    for lab, p in papers.items():
        for k, t in p["bib"].items():
            if k in SELF or k == SELF_SERIES:
                continue                      # unpublished parts of this corpus: cited as chapters
            if k not in bib:
                bib[k] = t
            elif bib[k] != t:
                if bib[k][:25].lower() == t[:25].lower():
                    if len(t) > len(bib[k]):
                        bib[k] = t
                    continue
                same = next((kk for kk, tt in bib.items() if tt == t), None)
                if same is None:
                    bib[f"{k}:{lab}"] = t
                    remap[(lab, k)] = f"{k}:{lab}"
                else:
                    remap[(lab, k)] = same
    parts_tex = []
    concord = []
    chno = 1   # introduction is chapter 1
    for part_title, chs in PLAN:
        parts_tex.append(f"\\part{{{pdfsafe(part_title)}}}")
        for folder, tex, lab in chs:
            chno += 1
            p = papers[lab]
            body = rewrite_cites(p["body"], {k: v for (l, k), v in remap.items() if l == lab})
            body = rewrite_paper_refs(body)
            body = re.sub(r"\\label\{", f"\\\\label{{{lab}:", body)
            body = re.sub(r"\\(ref|eqref|pageref|autoref)\{(?!chap:)", lambda m: "\\" + m.group(1) + "{" + lab + ":", body)
            body = re.sub(r"\\section\*\{Provenance\}", r"\\section*{Provenance of this chapter}", body)
            body = wrap_section_titles(body)
            # local macro overrides
            over = [text.replace("\\newcommand", "\\renewcommand", 1)
                    for name, text in p["macros"].items() if norm(text) != norm(macros[name])]
            # the table of contents carries the paper numeral and the short title
            short = re.sub(r"^Algorithmic non-commutative class field theory, ([IVX]+):\s*", "", p["title"])
            toc = (f"{lab}. " + short[0].upper() + short[1:]) if short != p["title"] else p["title"]
            fn = f"chapters/{lab}.tex"
            with open(os.path.join(outdir, fn), "w", encoding="utf-8") as f:
                f.write(f"% generated from {folder}/{tex}\n")
                f.write("\n".join(over) + ("\n" if over else ""))
                f.write(f"\\chapter[{pdfsafe(toc)}]{{{pdfsafe(p['title'])}}}\\label{{chap:{lab}}}\n")
                f.write(f"\\chaptermark{{{lab}}}\n")
                if p["abstract"]:
                    abst = re.sub(r"\\(ref|eqref|pageref|autoref)\{(?!chap:)", lambda m: "\\" + m.group(1) + "{" + lab + ":",
                                  rewrite_paper_refs(rewrite_cites(p["abstract"], {})))
                    f.write("\\begin{chapterabstract}\n" + abst + "\n\\end{chapterabstract}\n\n")
                f.write(body.strip() + "\n")
            parts_tex.append(f"\\input{{{fn[:-4]}}}")
            concord.append((lab, folder, chno, p["title"]))
    with open(os.path.join(outdir, "macros.tex"), "w", encoding="utf-8") as f:
        f.write("% union of the papers' macros; first definition wins, chapters override locally\n")
        for name, text in macros.items():
            f.write(text + "\n")
    with open(os.path.join(outdir, "bibliography.tex"), "w", encoding="utf-8") as f:
        f.write("\\begin{thebibliography}{999}\n")
        for k, t in bib.items():
            f.write(f"\\bibitem{{{k}}} {t}\n\n")
        f.write("\\end{thebibliography}\n")
    with open(os.path.join(outdir, "frontmatter", "concordance.tex"), "w", encoding="utf-8") as f:
        f.write("\\chapter*{Concordance}\n\\addcontentsline{toc}{chapter}{Concordance}\n"
                "The chapters are the papers of the series, re-ordered by subject. A reference of the form "
                "``[Ch.~$n$, Thm.~$a.b$]'' points to Theorem $n.a.b$ of this book; the papers' own "
                "section-wise numbering is preserved inside each chapter. References to the applied notes "
                "(Phys~1, Holo~1, Net~1, Ctrl~1) are by the labels below; the biological notes Bio~1--11 "
                "form Volume~II, \\emph{Algebraic Genomics}, and are not chapters of this volume.\n\n"
                "\\begin{longtable}{llr>{\\raggedright\\arraybackslash}p{8.2cm}}\n\\toprule\n"
                "paper & folder & chapter & title\\\\\n\\midrule\n\\endhead\n")
        for lab, folder, ch, title in concord:
            short = re.sub(r"^Algorithmic non-commutative class field theory, [IVX]+:\s*", "", title)
            f.write(f"{lab} & \\texttt{{{folder}}} & {ch} & {short}\\\\\n")
        f.write("\\bottomrule\n\\end{longtable}\n")
    with open(os.path.join(outdir, "main.tex"), "w", encoding="utf-8") as f:
        f.write(MAIN_TEMPLATE.replace("%%PARTS%%", "\n".join(parts_tex)))
    print(f"wrote {len(concord)} chapters, {len(macros)} macros, {len(bib)} references to {outdir}")


if __name__ == "__main__":
    papers = collections.OrderedDict()
    for folder, tex, lab in ORDER:
        papers[lab] = parse(os.path.join(locate(folder), tex))
    if sys.argv[1:] and sys.argv[1] == "analyse":
        analyse(papers)
    elif sys.argv[1:] and sys.argv[1] == "build":
        build(papers, sys.argv[2])
