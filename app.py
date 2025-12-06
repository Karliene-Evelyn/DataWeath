import os
import io
import tempfile
from datetime import datetime, date 

import requests
from requests.adapters import HTTPAdapter, Retry
import streamlit as st
import xarray as xr
import plotly.express as px
import plotly.graph_objects as go
import numpy as np

# Cartopy
try:
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.shapereader as shpreader
    HAS_CARTOPY = True
except Exception:
    HAS_CARTOPY = False
    shpreader = None

# Shapely
try:
    from shapely.geometry import Polygon, MultiPolygon, GeometryCollection, LineString, MultiLineString, box
    HAS_SHAPELY = True
except Exception:
    HAS_SHAPELY = False
    Polygon = MultiPolygon = GeometryCollection = LineString = MultiLineString = box = None

# ==== DESCRIÇÕES ====
from descricoes import descricao_isobaric, descricao_surface, descricao_meanSea

# =================== CONFIG PÁGINA ===================
st.set_page_config(page_title="DataWeath", layout="wide")

# ======= CSS =======
st.markdown("""
<style>
.block-container{ padding-top: 3rem !important; padding-bottom: 1rem !important; max-width: 1200px; }
h1.dataweath-title{ text-align:center; font-size:2.2rem; line-height:1.25; margin:0 0 .55rem 0; white-space:normal; word-break:break-word; overflow-wrap:anywhere; }
p.dataweath-subtitle{ text-align:center; color:#556; font-size:1.05rem; margin:.1rem 0 .9rem 0; }
.card{ background: var(--secondary-background-color); border:1px solid rgba(17,24,39,.06); border-radius:14px; padding:12px 14px; }
.chips{ display:flex; gap:.5rem; flex-wrap:wrap; margin:.2rem 0 .3rem 0 }
.chip{ background:#eef2ff; color:#3730a3; border-radius:999px; padding:.18rem .6rem; font-size:.8rem; border:1px solid rgba(55,48,163,.15)}
.btn-stack button{ width:100%; margin:.3rem 0; border-radius:12px !important; padding:.6rem 1rem !important; }
.stTabs [data-baseweb="tab-list"]{gap:6px; margin-top:.2rem; margin-bottom:.2rem;}
.stTabs [data-baseweb="tab"]{border-radius:12px; padding:.45rem .85rem}
.plotly-chart-container{ margin-top:.5rem; }
.small{ font-size:.92rem; color:#334155; }
</style>
""", unsafe_allow_html=True)

# ======= TÍTULO =======
st.markdown("""
<h1 class="dataweath-title">🌎 DataWeath - GFS Downloader e Viewer</h1>
<p class="dataweath-subtitle">Baixe, recorte e visualize dados GFS.</p>
""", unsafe_allow_html=True)

# =================== CONSTANTES ===================
NOAA_BASE = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
UCAR_BASE  = "https://osdf-data.gdex.ucar.edu/ncar/gdex/d084001"

NOAA_HINT  = "NOAA/NOMADS: rodadas recentes (últimas ~2 semanas)."
UCAR_HINT  = "UCAR/NCAR GDEX: histórico a partir de 15/01/2015."

OPEN_KW   = dict(engine="cfgrib", backend_kwargs={"indexpath": ""})
STEP_TYPES = ["instant", "avg", "accum"]
BORDERS_FALLBACK_URL = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_admin_0_boundary_lines_land.geojson"

# =================== HTTP SESSION ===================
@st.cache_resource
def get_http():
    s = requests.Session()
    s.headers.update({"User-Agent": "DataWeath/1.0"})
    retries = Retry(total=3, backoff_factor=0.6, status_forcelist=[429, 500, 502, 503, 504])
    s.mount("https://", HTTPAdapter(max_retries=retries))
    return s

def head_ok(url: str, timeout=(5, 10)):
    try:
        r = get_http().head(url, timeout=timeout, allow_redirects=True)
        return r.status_code == 200, r.headers
    except Exception:
        return False, {}

# =================== BUILD URLs ===================
def build_noaa_url_0p25(data_dt: datetime, hora: str, previsao: str) -> str:
    ymd = data_dt.strftime("%Y%m%d")
    ffff = previsao.lower().replace("f", "").zfill(3)
    return f"{NOAA_BASE}/gfs.{ymd}/{hora}/atmos/gfs.t{hora}z.pgrb2.0p25.f{ffff}"

def build_ucar_url_0p25(data_dt: datetime, hora: str, previsao: str) -> str:
    y = data_dt.strftime("%Y")
    ymd = data_dt.strftime("%Y%m%d")
    ffff = previsao.lower().replace("f", "").zfill(3)
    return f"{UCAR_BASE}/{y}/{ymd}/gfs.0p25.{ymd}{hora}.f{ffff}.grib2"

# =================== ESCOLHA AUTOMÁTICA PARA DOWNLOAD DE DADOS===================
def choose_best_source(data_dt: datetime, hora: str, previsao: str):
    noaa_url = build_noaa_url_0p25(data_dt, hora, previsao)
    ok, headers = head_ok(noaa_url)
    if ok:
        return "NOAA/NOMADS", noaa_url, headers
    ucar_url = build_ucar_url_0p25(data_dt, hora, previsao)
    ok2, headers2 = head_ok(ucar_url)
    if ok2:
        return "UCAR/NCAR GDEX", ucar_url, headers2
    return None, None, {}

# =================== SIDEBAR ===================
st.sidebar.header("📅 Parâmetros do GFS (0.25°)")
data_dt = st.sidebar.date_input(
    "Data da previsão",
    value=date(2015, 1, 15),
    min_value=date(2015, 1, 15),
    max_value=date.today(),
    help="Data da análise do modelo."
)
hora = st.sidebar.selectbox("Hora da análise (UTC)", ["00", "06", "12", "18"], help="Hora da análise de dados do GFS.")
previsao = st.sidebar.selectbox(
    "Horas à frente(Previsão)",
    ["f000", "f006", "f012", "f024", "f048", "f072", "f096", "f120", "f168", "f240"],
    help="Tempo de previsão a ser baixado."
)

st.sidebar.header("Recorte Espacial (Lat/Lon)")
colA, colB = st.sidebar.columns(2)
lat_min = colA.number_input("Lat mín (°)", value=-35.0, step=0.5)
lat_max = colB.number_input("Lat máx (°)", value=10.0, step=0.5)
colC, colD = st.sidebar.columns(2)
lon_min = colC.number_input("Lon mín (°)", value=-80.0, step=0.5)
lon_max = colD.number_input("Lon máx (°)", value=-34.0, step=0.5)

if lat_min >= lat_max or lon_min >= lon_max:
    st.sidebar.error("Verifique: mín < máx.")

with st.sidebar.expander("ℹ️ Ajuda", expanded=False):
    st.markdown(f"""
- **Fonte automática:** O sistema verifica automaticamente a origem dos dados:
                NOAA/NOMADS é a fonte prioritária com dados mais recentes.
                UCAR/NCAR GDEX é usada como alternativa caso a fonte principal esteja indisponível.
- **Disponibilidade:**
  - {NOAA_HINT}
  - {UCAR_HINT}
    """)

# ======= ESCOLHA DE FONTE =======
source_name, selected_url, selected_headers = choose_best_source(data_dt, hora, previsao)
st.markdown(f"""
<div class='card'>
  <div style="font-size:1.05rem"><b>🔗 Endpoint selecionado</b></div>
  <div class='chips'>
    <span class='chip'>Fonte: {source_name if source_name else "—"}</span>
    <span class='chip'>Resolução: 0.25°</span>
    <span class='chip'>Formato: GRIB2</span>
  </div>
  <div class="small">{selected_url if selected_url else "Nenhum endpoint disponível."}</div>
</div>
""", unsafe_allow_html=True)

# =================== DOWNLOAD ===================
def download_file(url: str):
    try:
        r = get_http().get(url, stream=True, timeout=(10, 120))
        if r.status_code >= 400:
            raise RuntimeError(f"HTTP {r.status_code}")
        total = int(r.headers.get("content-length", 0))
        pb = st.progress(0.0); done = 0
        fname = os.path.basename(url.split("?")[0]) or "gfs_0p25.grib2"
        temp_path = os.path.join(tempfile.gettempdir(), fname)
        with open(temp_path, "wb") as f:
            for ch in r.iter_content(1024 * 512):
                if not ch: continue
                f.write(ch); done += len(ch)
                if total: pb.progress(min(done / total, 1.0))
        st.success(f"✅ Baixado: {temp_path}")
        return temp_path, r.headers
    except Exception as e:
        st.error(f"Falha: {e}")
        return None, {}

# ====== Botões ======
st.markdown("<div class='btn-stack'>", unsafe_allow_html=True)
if st.button("🔍 Verificar arquivo"):
    if not selected_url:
        st.error("🚫 Nenhuma fonte válida. Tente outra data/hora.")
    else:
        ok, headers = head_ok(selected_url)
        st.session_state["url_ok"] = ok
        st.session_state["url_headers"] = headers
        if ok:
            size = headers.get("content-length")
            src = source_name or "desconhecida"
            info = f"Fonte: {src} • OK ✅"
            if size:
                info += f" • {int(size)/(1024*1024):.2f} MB"
            st.info(info)
        else:
            st.error("🚫 Link indisponível.")

can_dl = st.session_state.get("url_ok", False)
if st.button("📥 Baixar", disabled=not (can_dl and selected_url)):
    path_file, headers = download_file(selected_url)
    if path_file:
        st.session_state["path_file"] = path_file
        with open(path_file, "rb") as f:
            st.download_button("💾 Baixar GRIB2", f, file_name=os.path.basename(path_file))

st.markdown("</div>", unsafe_allow_html=True)

st.markdown("</div>", unsafe_allow_html=True)

# =================== FUNÇÕES DE DADOS ===================
@st.cache_data(show_spinner="🔍 Lendo variáveis disponíveis...", ttl=3600)
def listar_variaveis_e_niveis(path_file: str):
    variaveis = {}; niveis = {}
    # isobaric
    try:
        ds_iso = xr.open_dataset(path_file, filter_by_keys={"typeOfLevel": "isobaricInhPa"}, **OPEN_KW)
        variaveis["isobaricInhPa"] = [v for v in ds_iso.data_vars if ds_iso[v].ndim > 1]
        niveis["isobaricInhPa"] = list(ds_iso["isobaricInhPa"].values)
        ds_iso.close()
    except Exception:
        variaveis["isobaricInhPa"] = []; niveis["isobaricInhPa"] = []
    # surface
    surf_vars = set()
    for step in STEP_TYPES:
        try:
            ds_surface = xr.open_dataset(path_file, filter_by_keys={"typeOfLevel": "surface", "stepType": step}, **OPEN_KW)
            surf_vars |= {v for v in ds_surface.data_vars if ds_surface[v].ndim > 1}
            ds_surface.close()
        except Exception:
            pass
    variaveis["surface"] = sorted(surf_vars); niveis["surface"] = ["surface"]
    # meanSea
    try:
        ds_ms = xr.open_dataset(path_file, filter_by_keys={"typeOfLevel": "meanSea"}, **OPEN_KW)
        variaveis["meanSea"] = [v for v in ds_ms.data_vars if ds_ms[v].ndim > 1]
        niveis["meanSea"] = ["meanSea"]
        ds_ms.close()
    except Exception:
        variaveis["meanSea"] = []; niveis["meanSea"] = []
    return variaveis, niveis

def carregar_dataset(path_file: str, variavel: str, tipo_nivel: str, nivel: float | None = None) -> xr.DataArray:
    if tipo_nivel == "isobaricInhPa":
        ds = xr.open_dataset(path_file, filter_by_keys={"typeOfLevel": "isobaricInhPa"}, **OPEN_KW)
        try:
            if nivel not in ds["isobaricInhPa"].values:
                raise ValueError("Nível não disponível nesse arquivo.")
            da = ds[variavel].sel(isobaricInhPa=nivel)
            return da.load()
        finally:
            ds.close()
    elif tipo_nivel == "meanSea":
        ds = xr.open_dataset(path_file, filter_by_keys={"typeOfLevel": "meanSea"}, **OPEN_KW)
        try:
            if variavel not in ds.data_vars:
                raise ValueError(f"❌ Variável {variavel} não encontrada no nível meanSea.")
            return ds[variavel].load()
        finally:
            ds.close()
    else:  # surface
        for step in STEP_TYPES:
            ds = None
            try:
                ds = xr.open_dataset(path_file, filter_by_keys={"typeOfLevel": "surface", "stepType": step}, **OPEN_KW)
                if variavel in ds.data_vars:
                    return ds[variavel].load()
            except Exception:
                pass
            finally:
                try:
                    if ds is not None:
                        ds.close()
                except Exception:
                    pass
        raise ValueError(f"Variável {variavel} não encontrada em nenhum stepType do nível surface.")

# =================== RECORTE (BBOX) ===================
def normalize_lon_to_180(lon_da):
    return ((lon_da + 180) % 360) - 180

def subset_bbox(da: xr.DataArray, lat_min: float, lat_max: float, lon_min: float, lon_max: float) -> xr.DataArray:
    if da.coords["latitude"].values[0] > da.coords["latitude"].values[-1]:
        da = da.sortby("latitude")
    if float(da.coords["longitude"].max()) > 180:
        lon2 = normalize_lon_to_180(da.coords["longitude"])
        da = da.assign_coords(longitude=lon2).sortby("longitude")
    da_sub = da.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))
    return da_sub.load()

# =================== MAPA (Plotly + Fronteiras no bbox) ===================
def robust_limits(da, low=2, high=98):
    arr = np.asarray(da.values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None, None
    vmin = np.percentile(arr, low)
    vmax = np.percentile(arr, high)
    if np.isclose(vmin, vmax):
        return None, None
    return float(vmin), float(vmax)

def _extract_lines_from_geometry(geom):
    lines = []
    if HAS_SHAPELY:
        if isinstance(geom, Polygon):
            x, y = geom.exterior.coords.xy
            lines.append((np.asarray(x), np.asarray(y)))
        elif isinstance(geom, (MultiPolygon, GeometryCollection)):
            for g in geom.geoms:
                lines.extend(_extract_lines_from_geometry(g))
        elif isinstance(geom, LineString):
            x, y = geom.coords.xy
            lines.append((np.asarray(x), np.asarray(y)))
        elif isinstance(geom, MultiLineString):
            for g in geom.geoms:
                x, y = g.coords.xy
                lines.append((np.asarray(x), np.asarray(y)))
    return lines

def _clip_xy_to_bbox(x, y, bbox_poly):
    if not HAS_SHAPELY:
        return [(x, y)]
    try:
        line = LineString(np.column_stack([x, y]))
        inter = line.intersection(bbox_poly)
        clipped = []
        if inter.is_empty:
            return []
        if isinstance(inter, (LineString,)):
            xs, ys = inter.coords.xy
            clipped.append((np.asarray(xs), np.asarray(ys)))
        elif isinstance(inter, MultiLineString):
            for seg in inter.geoms:
                xs, ys = seg.coords.xy
                clipped.append((np.asarray(xs), np.asarray(ys)))
        return clipped
    except Exception:
        return [(x, y)]

@st.cache_data(show_spinner=False)
def _fetch_borders_geojson():
    try:
        r = requests.get(BORDERS_FALLBACK_URL, timeout=20)
        r.raise_for_status()
        return r.json()
    except Exception:
        return None

def _geojson_lines_clipped(geojson_obj, lon_min, lon_max, lat_min, lat_max):
    traces = []
    bbox_poly = box(lon_min, lat_min, lon_max, lat_max) if HAS_SHAPELY else None
    try:
        feats = geojson_obj.get("features", [])
        for f in feats:
            geom = f.get("geometry", {})
            gtype = geom.get("type")
            coords = geom.get("coordinates")
            if gtype == "LineString":
                xs, ys = zip(*coords)
                xs = np.asarray(xs); ys = np.asarray(ys)
                if (xs.min() > lon_max) or (xs.max() < lon_min) or (ys.min() > lat_max) or (ys.max() < lat_min):
                    continue
                if HAS_SHAPELY:
                    for xs_c, ys_c in _clip_xy_to_bbox(xs, ys, bbox_poly):
                        traces.append((xs_c, ys_c))
                else:
                    traces.append((xs, ys))
            elif gtype == "MultiLineString":
                for seg in coords:
                    xs, ys = zip(*seg)
                    xs = np.asarray(xs); ys = np.asarray(ys)
                    if (xs.min() > lon_max) or (xs.max() < lon_min) or (ys.min() > lat_max) or (ys.max() < lat_min):
                        continue
                    if HAS_SHAPELY:
                        for xs_c, ys_c in _clip_xy_to_bbox(xs, ys, bbox_poly):
                            traces.append((xs_c, ys_c))
                    else:
                        traces.append((xs, ys))
    except Exception:
        pass
    return traces

def get_border_traces_for_bbox(lon_min, lon_max, lat_min, lat_max,
                               linecolor_coast="rgba(0,0,0,0.45)", lw_coast=0.8,
                               linecolor_bnd="rgba(0,0,0,0.35)", lw_bnd=0.6):
    traces = []
    if HAS_CARTOPY and shpreader is not None:
        try:
            bbox_poly = box(lon_min, lat_min, lon_max, lat_max) if HAS_SHAPELY else None
            coast_path = shpreader.natural_earth(resolution="50m", category="physical", name="coastline")
            for rec in shpreader.Reader(coast_path).records():
                geom = rec.geometry
                geoms = [geom]
                if HAS_SHAPELY and bbox_poly is not None:
                    inter = geom.intersection(bbox_poly)
                    geoms = [inter] if not isinstance(inter, (MultiPolygon, GeometryCollection)) else list(inter.geoms)
                for g in geoms:
                    for x, y in _extract_lines_from_geometry(g):
                        traces.append(go.Scatter(x=x, y=y, mode="lines",
                                                 line=dict(color=linecolor_coast, width=lw_coast),
                                                 hoverinfo="skip", showlegend=False))
            bnd_path = shpreader.natural_earth(resolution="50m", category="cultural", name="admin_0_boundary_lines_land")
            for rec in shpreader.Reader(bnd_path).records():
                geom = rec.geometry
                geoms = [geom]
                if HAS_SHAPELY and bbox_poly is not None:
                    inter = geom.intersection(bbox_poly)
                    geoms = [inter] if not isinstance(inter, (MultiPolygon, GeometryCollection)) else list(inter.geoms)
                for g in geoms:
                    for x, y in _extract_lines_from_geometry(g):
                        traces.append(go.Scatter(x=x, y=y, mode="lines",
                                                 line=dict(color=linecolor_bnd, width=lw_bnd),
                                                 hoverinfo="skip", showlegend=False))
            return traces
        except Exception:
            traces = []
    # Fallback GeoJSON
    geojson_obj = _fetch_borders_geojson()
    if geojson_obj:
        for x, y in _geojson_lines_clipped(geojson_obj, lon_min, lon_max, lat_min, lat_max):
            traces.append(go.Scatter(x=x, y=y, mode="lines",
                                     line=dict(color=linecolor_bnd, width=lw_bnd),
                                     hoverinfo="skip", showlegend=False))
    return traces

def make_plotly_imshow_clean(dados, title, colorscale="Viridis", add_borders=True, add_light_contours=False, height=520):
    lat = dados.coords["latitude"].values
    lon = dados.coords["longitude"].values
    z = dados.values
    if lat[0] > lat[-1]:
        lat = lat[::-1]; z = z[::-1, :]
    vmin, vmax = robust_limits(dados)
    fig = px.imshow(
        z, x=lon, y=lat, origin="lower",
        labels={"x": "Longitude", "y": "Latitude", "color": dados.attrs.get("long_name", dados.name)},
        aspect="auto", color_continuous_scale=colorscale, zmin=vmin, zmax=vmax
    )
    fig.update_layout(
        title=title, height=height, margin=dict(l=10, r=10, t=28, b=6),
        coloraxis_colorbar=dict(title=dados.attrs.get("units", ""), thickness=14, len=0.8),
    )
    fig.update_xaxes(showgrid=False, ticks="outside", ticklen=4, tickwidth=1, mirror=True)
    fig.update_yaxes(showgrid=False, ticks="outside", ticklen=4, tickwidth=1, mirror=True)

    if add_light_contours:
        fig.add_contour(z=z, x=lon, y=lat, contours=dict(showlabels=False, coloring="none"),
                        line_width=0.6, line_color="rgba(0,0,0,0.25)", showscale=False)

    if add_borders:
        lon_min_v, lon_max_v = float(lon.min()), float(lon.max())
        lat_min_v, lat_max_v = float(lat.min()), float(lat.max())
        for tr in get_border_traces_for_bbox(lon_min_v, lon_max_v, lat_min_v, lat_max_v):
            fig.add_trace(tr)
    return fig

# =================== CARTOPY (HQ) ===================
def make_cartopy_fig(dados, title, draw_contours=False):
    if not HAS_CARTOPY:
        st.warning("Cartopy não instalado — usando apenas Plotly.")
        return None
    da = dados
    if da.coords["latitude"].values[0] > da.coords["latitude"].values[-1]:
        da = da.sortby("latitude")
    if float(da.coords["longitude"].max()) > 180:
        lon2 = ((da.coords["longitude"] + 180) % 360) - 180
        da = da.assign_coords(longitude=lon2).sortby("longitude")

    lats = da.latitude.values
    lons = da.longitude.values
    Z = da.values
    vmin, vmax = robust_limits(da)

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10.5, 6.8), dpi=220)
    ax = plt.axes(projection=proj)

    pm = ax.pcolormesh(lons, lats, Z, transform=proj, cmap="viridis",
                       shading="auto", vmin=vmin, vmax=vmax)
    cb = plt.colorbar(pm, ax=ax, shrink=0.90, pad=0.03)
    cb.set_label(da.attrs.get("units", ""))

    if draw_contours:
        try:
            ax.contour(lons, lats, Z, levels=8, colors="k", linewidths=0.35, alpha=0.35, transform=proj)
        except Exception:
            pass

    ax.coastlines(resolution="50m", linewidth=0.7)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), edgecolor="k", facecolor="none", linewidth=0.3)

    gl = ax.gridlines(draw_labels=True, linestyle="--", color="gray", alpha=0.35, linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    ax.set_extent([float(da.longitude.min()), float(da.longitude.max()),
                   float(da.latitude.min()), float(da.latitude.max())], crs=proj)

    ax.set_title(title, fontsize=12, pad=6)
    fig.tight_layout()
    return fig

# =================== TABS (Visualização • Biblioteca) ===================
tab_viz, tab_lib = st.tabs(["📊 Visualizar", "📚 Biblioteca"])

with tab_viz:
    if "path_file" in st.session_state:
        variaveis_dict, niveis_dict = listar_variaveis_e_niveis(st.session_state["path_file"])

        # ---- SEM AUTOSELEÇÃO ----
        tipo_nivel_opcoes = {
            "isobaricInhPa": "Camadas de Pressão (hPa)",
            "surface": "Superfície",
            "meanSea": "Nível Médio do Mar (MSLP)",
        }
        tipo_nivel_legivel = st.selectbox("Tipo de nível", ["Selecione..." ] + list(tipo_nivel_opcoes.values()), index=0, help="Escolha o nível vertical para listar as variáveis disponíveis.")

        tipo_nivel_escolhido = None
        for key, label in tipo_nivel_opcoes.items():
            if tipo_nivel_legivel == label:
                tipo_nivel_escolhido = key

        if tipo_nivel_escolhido:
            variaveis = variaveis_dict.get(tipo_nivel_escolhido, [])
            if variaveis:
                var_escolhida = st.selectbox("Variável", ["Selecione..."] + variaveis, index=0, help="Selecione a variável meteorológica que deseja visualizar.")
                if var_escolhida != "Selecione...":
                    # descrição
                    if tipo_nivel_escolhido == "isobaricInhPa":
                        descricao_ativa = descricao_isobaric
                    elif tipo_nivel_escolhido == "surface":
                        descricao_ativa = descricao_surface
                    elif tipo_nivel_escolhido == "meanSea":
                        descricao_ativa = descricao_meanSea
                    else:
                        descricao_ativa = {}
                    if var_escolhida in descricao_ativa:
                        st.caption(descricao_ativa[var_escolhida])

                    # nível (isobárico)
                    if tipo_nivel_escolhido == "isobaricInhPa":
                        niveis_legiveis = [str(n) for n in niveis_dict.get("isobaricInhPa", [])]
                        nivel_pick = st.selectbox("Nível (hPa)", ["Selecione..."] + niveis_legiveis, index=0, help="Pressão em hPa (ex.: 850, 500).")
                        if nivel_pick != "Selecione...":
                            dados = carregar_dataset(st.session_state["path_file"], var_escolhida, "isobaricInhPa", float(nivel_pick))
                        else:
                            dados = None
                    else:
                        dados = carregar_dataset(st.session_state["path_file"], var_escolhida, tipo_nivel_escolhido)

                    # recorte + plot
                    if dados is not None:
                        try:
                            dados = subset_bbox(dados, lat_min, lat_max, lon_min, lon_max)
                        except Exception as e:
                            st.warning(f"Não foi possível recortar pelo bbox: {e}")

                        if dados.sizes.get("latitude", 0) < 2 or dados.sizes.get("longitude", 0) < 2:
                            st.warning("⚠️ Recorte muito pequeno para visualização (menos de 2x2 pontos).")
                        else:
                            titulo = f"{dados.name} ({dados.attrs.get('units','')}) • {data_dt.strftime('%Y-%m-%d')} {hora}Z • {previsao.upper()}"

                            tab1, tab2 = st.tabs(["⚡ Plotly (bbox)", "🖼️ Cartopy (HQ)"])

                            with tab1:
                                fig = make_plotly_imshow_clean(
                                    dados, titulo, colorscale="Viridis",
                                    add_borders=True, add_light_contours=False, height=520
                                )
                                st.markdown('<div class="plotly-chart-container">', unsafe_allow_html=True)
                                st.plotly_chart(fig, use_container_width=True)
                                st.markdown('</div>', unsafe_allow_html=True)

                            with tab2:
                                if HAS_CARTOPY:
                                    draw_ct = st.checkbox("Contornos leves no Cartopy", value=False, help="Desenhar isolinhas suaves por cima do sombreado.")
                                    figc = make_cartopy_fig(dados, titulo, draw_contours=draw_ct)
                                    if figc is not None:
                                        buf = io.BytesIO()
                                        figc.savefig(buf, format="png", dpi=220, bbox_inches="tight")
                                        buf.seek(0)
                                        st.download_button("📸 Baixar PNG (HQ)", buf, file_name="mapa_hq.png", mime="image/png")
                                        st.pyplot(figc, clear_figure=True)
                                        plt.close(figc)
                                else:
                                    st.info("Cartopy não instalado — a aba HQ fica desativada no seu ambiente.")
            else:
                st.warning("Nenhuma variável disponível para esse tipo de nível.")
        else:
            st.info("Escolha um tipo de nível para continuar.")
    else:
        st.info("Verifique e baixe o arquivo acima para habilitar a visualização.")

with tab_lib:
    st.subheader("Biblioteca de Variáveis")
    nivel_opcoes = {
        "isobaricInhPa": "Camadas de Pressão (hPa)",
        "surface": "Superfície",
        "meanSea": "Nível Médio do Mar (MSLP)",
    }
    colN1, colN2 = st.columns([1, 1.2])
    nivel_legivel = colN1.selectbox("Tipo de nível:", list(nivel_opcoes.values()), help="Filtre o dicionário pelo tipo de nível.")
    query = colN2.text_input("🔎 Buscar por nome/descrição", "", help="Procure pela sigla da variável ou por termos da descrição.")

    nivel_escolhido = None
    for key, label in nivel_opcoes.items():
        if nivel_legivel == label:
            nivel_escolhido = key

    if nivel_escolhido:
        if nivel_escolhido == "isobaricInhPa":
            dicionario = descricao_isobaric
        elif nivel_escolhido == "surface":
            dicionario = descricao_surface
        elif nivel_escolhido == "meanSea":
            dicionario = descricao_meanSea
        else:
            dicionario = {}

        if query:
            q = query.lower().strip()
            dicionario = {k: v for k, v in dicionario.items() if q in k.lower() or q in v.lower()}

        if dicionario:
            st.markdown("<div class='card'>", unsafe_allow_html=True)
            dicionario_ordenado = dict(sorted(dicionario.items(), key=lambda item: item[0].lower()))
            for var, desc in dicionario_ordenado.items():
                st.markdown(f"**`{var}`** — {desc}")
            st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.info("Nada encontrado para esse filtro.")