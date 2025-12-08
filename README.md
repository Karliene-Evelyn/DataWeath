# 🌎 DataWeath — GFS Downloader e Viewer

O **DataWeath** é uma aplicação em Streamlit para baixar, recortar e visualizar dados do modelo meteorológico **GFS 0.25°**.  
Oferece dois modos de visualização:

- **Plotly (interativo)** — funciona em qualquer ambiente.
- **Cartopy (alta resolução)** — funciona apenas **localmente**, pois depende de bibliotecas de sistema que não existem no Streamlit Cloud.

---

## 1. Instalação

### 1.1. Clonar o repositório

```bash
git clone https://github.com/Karliene-Evelyn/DataWeath.git
cd DataWeath
```
1.2. Criar ambiente virtual (opcional, mas recomendado)
```bash
python -m venv .venv
```
# Windows
```bash
.venv\Scripts\activate
```
# Linux/macOS
```bash
source .venv/bin/activate
```
1.3. Instalar bibliotecas principais
```
pip install streamlit requests xarray plotly numpy shapely matplotlib
pip install cfgrib
```
Obs.: O cfgrib precisa do ecCodes instalado na máquina.
Se preferir instalar tudo automaticamente (incluindo Cartopy e ecCodes), use Conda:
```
conda install -c conda-forge xarray cfgrib eccodes cartopy shapely matplotlib plotly requests streamlit
```
✔️ 2. Executar a aplicação
No diretório do projeto:
```bash
streamlit run app.py
```
A aplicação abrirá no navegador em:
```
http://localhost:8501
```
✔️ 3. Estrutura básica do projeto
```bash
DataWeath/
├─ app.py                # Aplicação principal
├─ descricoes.py         # Dicionários com descrições de variáveis
└─ README.md
