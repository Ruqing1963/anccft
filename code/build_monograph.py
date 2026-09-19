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

# (part title, [(folder, tex, roman/label)])  -- the monograph plan
PLAN = [
    ("The bridge over the tree: two channels and one obstruction", [
        ("Paper 1", "anccft.tex", "I"), ("Paper 2", "anccft2.tex", "II"),
        ("Paper 6", "anccft6.tex", "VI"), ("Paper 3", "anccft3.tex", "III")]),
    ("Towers: commutative Iwasawa theory and its failure", [
        ("Paper 4", "anccft4.tex", "IV"), ("Paper 9", "anccft9.tex", "IX"),
        ("Paper 14", "anccft14.tex", "XIV"), ("Paper 5", "anccft5.tex", "V"),
        ("Paper 7", "anccft7.tex", "VII"), ("Paper 13", "anccft13.tex", "XIII")]),
    ("Operator algebras of the tower", [
        ("Paper 8", "anccft8.tex", "VIII"), ("Paper 16", "anccft16.tex", "XVI")]),
    ("Rank two: $\\widetilde A_2$ buildings", [
        ("Paper 10", "anccft10.tex", "X"), ("Paper 11", "anccft11.tex", "XI"),
        ("Paper 12", "anccft12.tex", "XII"), ("Paper 15", "anccft15.tex", "XV")]),
    ("Applications outside number theory", [
        ("Phys 1", "bloch_mass.tex", "Phys1"), ("Holo 1", "boundary_blind.tex", "Holo1"),
        ("Net 1", "tower_lambda.tex", "Net1"), ("Ctrl 1", "hidden_eigs.tex", "Ctrl1"),
        ("Bio 1", "dbg_sandpile.tex", "Bio1"), ("Bio 2", "rotor_assembly.tex", "Bio2"),
        ("Bio 3", "repeat_splitting.tex", "Bio3"), ("Bio 4", "multicopy.tex", "Bio4")]),
]
ORDER = [(f, t, l) for _, chs in PLAN for f, t, l in chs]

MACRO_RE = re.compile(r"\\(?:re)?newcommand\*?\{?\\(\w+)\}?(\[\d\])?(\[[^\]]*\])?\{")
MATHOP_RE = re.compile(r"\\DeclareMathOperator\*?\{\\(\w+)\}\{([^}]*)\}")
THM_RE = re.compile(r"\\newtheorem\*?\{(\w+)\}(?:\[(\w+)\])?\{([^}]*)\}(?:\[(\w+)\])?")
ROMAN = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"}


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


def rewrite_cites(body, mapping):
    def rep(m):
        keys = [mapping.get(k.strip(), k.strip()) for k in m.group(2).split(",")]
        return m.group(1) + "{" + ",".join(keys) + "}"
    return re.sub(r"(\\cite[tp]?\*?(?:\[[^\]]*\])?)\{([^}]*)\}", rep, body)


def rewrite_paper_refs(body):
    """[V, Thm.~3.1] -> [Ch.~\ref{chap:V}, Thm.~3.1];  [V] -> [Ch.~\ref{chap:V}] when clearly a paper."""
    def rep(m):
        r = m.group(1)
        return f"[Ch.~\\ref{{chap:{r}}}, " if r in ROMAN else m.group(0)
    body = re.sub(r"\[([IVX]{1,5}),\s*(?=(?:Thm|Prop|Lem|Cor|Def|Rem|Prob|Tab|Eq|Conv|Constr|Conj|Alg|Ex|Sec|Fig|Hyp|eq|\\S|\\eqref|\\S|\(|\u00a7|Section|Table|Remark|Theorem))", rep, body)
    body = re.sub(r"(?<![\w$^_{\\])\[([IVX]{1,5})\](?=[\s.,;:)])",
                  lambda m: f"[Ch.~\\ref{{chap:{m.group(1)}}}]" if m.group(1) in ROMAN else m.group(0), body)
    # ranges and lists like [I]--[IV] or [X]--[XII]
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
\date{September 2026}

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
        parts_tex.append(f"\\part{{{part_title}}}")
        for folder, tex, lab in chs:
            chno += 1
            p = papers[lab]
            body = rewrite_cites(p["body"], {k: v for (l, k), v in remap.items() if l == lab})
            body = rewrite_paper_refs(body)
            body = re.sub(r"\\label\{", f"\\\\label{{{lab}:", body)
            body = re.sub(r"\\(ref|eqref|pageref|autoref)\{(?!chap:)", lambda m: "\\" + m.group(1) + "{" + lab + ":", body)
            body = re.sub(r"\\section\*\{Provenance\}", r"\\section*{Provenance of this chapter}", body)
            # local macro overrides
            over = [text.replace("\\newcommand", "\\renewcommand", 1)
                    for name, text in p["macros"].items() if norm(text) != norm(macros[name])]
            fn = f"chapters/{lab}.tex"
            with open(os.path.join(outdir, fn), "w", encoding="utf-8") as f:
                f.write(f"% generated from {folder}/{tex}\n")
                f.write("\n".join(over) + ("\n" if over else ""))
                f.write(f"\\chapter{{{p['title']}}}\\label{{chap:{lab}}}\n")
                f.write(f"\\chaptermark{{{lab}}}\n")
                if p["abstract"]:
                    abst = re.sub(r"\\(ref|eqref|pageref|autoref)\{(?!chap:)", lambda m: "\\" + m.group(1) + "{" + lab + ":",
                                  rewrite_paper_refs(p["abstract"]))
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
                "(Phys~1, Bio~1--4, \\dots) are by the labels below.\n\n"
                "\\begin{longtable}{llrp{8.2cm}}\n\\toprule\npaper & folder & chapter & title\\\\\n\\midrule\n\\endhead\n")
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
        papers[lab] = parse(os.path.join(ROOT, folder, tex))
    if sys.argv[1:] and sys.argv[1] == "analyse":
        analyse(papers)
    elif sys.argv[1:] and sys.argv[1] == "build":
        build(papers, sys.argv[2])
