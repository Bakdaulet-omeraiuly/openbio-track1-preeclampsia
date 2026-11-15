"""Build the narrated, animated presentation video from narration.json and the repo figures.

Narration: macOS `say` (offline text-to-speech). Animation: slow zoom + fades per slide (ffmpeg).
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
BG, INK, MUTED, BLUE = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6"
W, H = 1920, 1080
slides = json.load(open(HERE / "narration.json"))


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def render(s, path):
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=BG)
    fig.text(0.05, 0.92, s["title"], fontsize=40 if len(s["title"]) <= 50 else 34, color=INK, weight="bold", va="top")
    fig.add_artist(plt.Line2D([0.05, 0.25], [0.855, 0.855], color=BLUE, lw=6, transform=fig.transFigure))
    if s.get("sub"):
        fig.text(0.05, 0.825, s["sub"], fontsize=22, color=MUTED, va="top")
    if s.get("fig"):
        img = mpimg.imread(FIG / s["fig"])
        ih, iw = img.shape[:2]
        box_w, box_h = 0.9, 0.68
        scale = min(box_w * W / iw, box_h * H / ih)
        w, h = iw * scale / W, ih * scale / H
        ax = fig.add_axes([0.5 - w / 2, 0.06 + (box_h - h) / 2, w, h]); ax.imshow(img); ax.axis("off")
    if s.get("bullets"):
        y = 0.76
        for b in s["bullets"]:
            for k, line in enumerate(textwrap.wrap(b, 95)):
                fig.text(0.07, y, ("•  " if k == 0 else "    ") + line, fontsize=24, color=INK, va="top"); y -= 0.055
            y -= 0.025
    if s["id"] == "s01_title":
        fig.text(0.05, 0.55, "Placental gene activity, 10 patient groups, single cells, intact tissue and knockouts", fontsize=26, color=INK)
        fig.text(0.05, 0.47, "Every test pre-registered; positives and negatives reported", fontsize=24, color=MUTED)
    fig.savefig(path, facecolor=BG); plt.close(fig)


def dur(audio):
    r = subprocess.run([FF, "-i", str(audio)], capture_output=True, text=True).stderr
    t = r.split("Duration: ")[1].split(",")[0]; h, m, sec = t.split(":")
    return int(h) * 3600 + int(m) * 60 + float(sec)


clips = []
for s in slides:
    png, aiff, mp4 = OUT / f"{s['id']}.png", OUT / f"{s['id']}.aiff", OUT / f"{s['id']}.mp4"
    render(s, png)
    run(["say", "-v", "Samantha", "-r", "165", "-o", str(aiff), s["text"]])
    d = dur(aiff) + 1.2
    frames = int(d * 30)
    vf = (f"scale=3840:-1,zoompan=z='min(zoom+0.00025,1.06)':x='iw/2-(iw/zoom/2)':y='ih/2-(ih/zoom/2)':d={frames}:s={W}x{H}:fps=30,"
          f"fade=t=in:st=0:d=0.5,fade=t=out:st={d - 0.5:.2f}:d=0.5,format=yuv420p")
    run([FF, "-y", "-loop", "1", "-i", str(png), "-i", str(aiff), "-filter_complex",
         f"[0:v]{vf}[v];[1:a]adelay=400|400,loudnorm=I=-18:TP=-1.5:LRA=11,aresample=48000,apad[a]", "-map", "[v]", "-map", "[a]", "-t", f"{d:.2f}",
         "-c:v", "libx264", "-preset", "medium", "-crf", "20", "-c:a", "aac", "-b:a", "160k", "-ar", "48000", "-ac", "2", str(mp4)])
    clips.append(mp4); print(f"{s['id']}: {d:.1f}s")


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
