"""Build the narrated, animated presentation video from narration.json and the repo figures.

Narration: macOS `say` or edge-tts (VOICE=edge:en-US-AvaNeural; VOICE env var, e.g. 'Ava (Premium)'). Animation: staged reveals with crossfades (ffmpeg xfade).
Optional on-camera clips by the author: put intro.mov / outro.mov in this folder; they are joined at start/end.
Usage: python build_video.py FFMPEG_EXE
"""
import json, subprocess, sys, textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

FF = sys.argv[1]
HERE = Path(__file__).resolve().parent
FIG = HERE.parent / "figures"
OUT = HERE / "build"; OUT.mkdir(exist_ok=True)
import os
BG, INK, MUTED, BLUE = "#0e1f1a", "#f2f1ec", "#a9b5ae", "#f0a35e"
CARD = "#fcfcfb"
VOICE = os.environ.get("VOICE", "Samantha")
W, H = 1920, 1080
slides = json.load(open(HERE / "narration.json"))


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _base(s):
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=BG)
    fig.text(0.05, 0.92, s["title"], fontsize=44 if len(s["title"]) <= 46 else 36, color=INK, weight="bold", va="top")
    fig.add_artist(plt.Line2D([0.05, 0.17], [0.845, 0.845], color=BLUE, lw=7, transform=fig.transFigure, solid_capstyle="round"))
    if s.get("sub"):
        fig.text(0.05, 0.815, s["sub"], fontsize=22, color=MUTED, va="top")
    if s["id"] == "s01_title":
        fig.text(0.05, 0.60, "Placental gene activity · 10 patient groups · single cells · intact tissue · knockouts", fontsize=26, color=INK)
        fig.text(0.05, 0.52, "Every test pre-registered; positives and negatives reported", fontsize=24, color=BLUE)
    return fig


def _add_fig(fig, name):
    img = mpimg.imread(FIG / name)
    ih, iw = img.shape[:2]
    box_w, box_h = 0.86, 0.66
    scale = min(box_w * W / iw, box_h * H / ih)
    w, h = iw * scale / W, ih * scale / H
    x0, y0 = 0.5 - w / 2, 0.06 + (box_h - h) / 2
    from matplotlib.patches import FancyBboxPatch
    fig.patches.append(FancyBboxPatch((x0 - 0.012, y0 - 0.02), w + 0.024, h + 0.04, boxstyle="round,pad=0,rounding_size=0.015",
                                      transform=fig.transFigure, color=CARD, zorder=-5))
    ax = fig.add_axes([x0, y0, w, h], zorder=5); ax.imshow(img); ax.axis("off")


def render_states(s):
    """Cumulative frames: title -> (+figure) -> +bullet 1 -> +bullet 2 ... Returns list of PNG paths."""
    steps = [("base", None)]
    if s.get("fig"):
        steps.append(("fig", s["fig"]))
    for k in range(len(s.get("bullets", []))):
        steps.append(("bullet", k))
    paths = []
    for n in range(1, len(steps) + 1):
        fig = _base(s)
        y = 0.74
        nb = len(s.get("bullets", []))
        fs, wrap, gap = (25, 92, 0.056) if nb <= 5 else (21, 110, 0.047)
        for kind, val in steps[:n]:
            if kind == "fig":
                _add_fig(fig, val)
            elif kind == "bullet":
                b = s["bullets"][val]
                for k, line in enumerate(textwrap.wrap(b, wrap)):
                    fig.text(0.07, y, ("●  " if k == 0 else "     ") + line, fontsize=fs, color=INK, va="top"); y -= gap
                y -= gap / 2
        if kind == "bullet":
            pass
        path = OUT / f"{s['id']}_state{n}.png"; fig.savefig(path, facecolor=BG); plt.close(fig); paths.append(path)
    return paths


def dur(audio):
    r = subprocess.run([FF, "-i", str(audio)], capture_output=True, text=True).stderr
    t = r.split("Duration: ")[1].split(",")[0]; h, m, sec = t.split(":")
    return int(h) * 3600 + int(m) * 60 + float(sec)


clips = []
for s in slides:
    states = render_states(s)
    aiff, mp4 = OUT / f"{s['id']}.aiff", OUT / f"{s['id']}.mp4"
    if VOICE.startswith("edge:"):  # Microsoft neural voice via edge-tts (sends only the narration text)
        mp3 = aiff.with_suffix(".mp3")
        run([sys.executable, "-m", "edge_tts", "--voice", VOICE[5:], "--rate", "+12%", "--text", s["text"].replace("s F L T one", "soluble F-L-T-1").replace("F L T one", "F-L-T-1"), "--write-media", str(mp3)])
        aiff = mp3
    else:
        run(["say", "-v", VOICE, "-r", "170", "-o", str(aiff), s["text"]])
    d = dur(aiff) + 1.0
    n = len(states)
    # reveal times: first state at 0, the rest spread over the first 75% of the narration
    starts = [0.0] + [0.6 + (d * 0.75 - 0.6) * k / max(n - 1, 1) for k in range(1, n)]
    segs = [(starts[k + 1] if k + 1 < n else d) - starts[k] for k in range(n)]
    XF = 0.5
    cmd = [FF, "-y"]
    for k, (p, sd) in enumerate(zip(states, segs)):
        cmd += ["-loop", "1", "-t", f"{sd + (XF if k + 1 < n else 0):.3f}", "-i", str(p)]
    cmd += ["-i", str(aiff)]
    chain, last, off = [], "0:v", 0.0
    for k in range(1, n):
        off += segs[k - 1]
        chain.append(f"[{last}][{k}:v]xfade=transition=fade:duration={XF}:offset={off - XF + XF:.3f}[x{k}]" if False else
                     f"[{'0:v' if k == 1 else f'x{k-1}'}][{k}:v]xfade=transition=fade:duration={XF}:offset={off:.3f}[x{k}]")
    vlast = f"x{n-1}" if n > 1 else "0:v"
    fc = ";".join(chain + [f"[{vlast}]fps=30,format=yuv420p,fade=t=in:st=0:d=0.4,fade=t=out:st={d - 0.4:.2f}:d=0.4[v]",
                           f"[{n}:a]adelay=300|300,loudnorm=I=-18:TP=-1.5:LRA=11,aresample=48000,apad[a]"])
    cmd += ["-filter_complex", fc, "-map", "[v]", "-map", "[a]", "-t", f"{d:.2f}", "-c:v", "libx264", "-preset", "medium",
            "-crf", "19", "-c:a", "aac", "-b:a", "160k", "-ar", "48000", "-ac", "2", str(mp4)]
    run(cmd)
    clips.append(mp4); print(f"{s['id']}: {d:.1f}s, {n} states")


def norm(src, dst):
    """Author clip: HDR (HLG) phone video -> SDR; vertical video centred on a blurred copy of itself; first audio track."""
    tone = ("zscale=tin=arib-std-b67:min=bt2020nc:pin=bt2020:t=linear:npl=100,format=gbrpf32le,"
            "zscale=p=bt709,tonemap=hable:desat=0,zscale=t=bt709:m=bt709:r=tv,format=yuv420p")
    vf = (f"[0:v:0]{tone},fps=30,split[a][b];"
          f"[a]scale={W}:{H}:force_original_aspect_ratio=increase,crop={W}:{H},gblur=sigma=30,eq=brightness=-0.08[bg];"
          f"[b]scale=-2:{H}[fg];[bg][fg]overlay=(W-w)/2:0,format=yuv420p,fade=t=in:st=0:d=0.4[v]")
    run([FF, "-y", "-i", str(src), "-filter_complex", vf, "-map", "[v]", "-map", "0:a:0", "-af", "loudnorm=I=-18:TP=-1.5:LRA=11,aresample=48000",
         "-c:v", "libx264", "-crf", "20", "-preset", "medium", "-c:a", "aac", "-b:a", "160k", "-ar", "48000", "-ac", "2", str(dst)])
    return dst


seq = []
for name in ("intro",):
    for ext in (".mov", ".mp4", ".m4v"):
        if (HERE / f"{name}{ext}").exists():
            seq.append(norm(HERE / f"{name}{ext}", OUT / f"{name}_norm.mp4")); break
seq += clips
for name in ("outro",):
    for ext in (".mov", ".mp4", ".m4v"):
        if (HERE / f"{name}{ext}").exists():
            seq.append(norm(HERE / f"{name}{ext}", OUT / f"{name}_norm.mp4")); break
(OUT / "list.txt").write_text("".join(f"file '{p}'\n" for p in seq))
final = HERE / "presentation.mp4"
run([FF, "-y", "-f", "concat", "-safe", "0", "-i", str(OUT / "list.txt"), "-c:v", "libx264", "-crf", "20", "-c:a", "aac", "-b:a", "160k", str(final)])
print("total", round(dur(final), 1), "s ->", final, "| intro/outro included:", [p.name for p in seq if "norm" in p.name])
