# get_peptide_card.py
import sys
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

def get_peptide_card_info(peptide_id):
    url = f"https://dbaasp.org/peptide-card?id=DBAASPS_{peptide_id}"

    options = Options()
    # Modo headless moderno
    options.add_argument("--headless=new")
    # Para Windows: desactivar GPU / composición directa y 3D
    options.add_argument("--disable-gpu")
    options.add_argument("--disable-software-rasterizer")
    options.add_argument("--disable-3d-apis")
    options.add_argument("--disable-direct-composition")
    # Forzar backend de GL por software (evita errores de DirectComposition)
    options.add_argument("--use-angle=swiftshader")
    options.add_argument("--use-gl=swiftshader")

    # Estabilidad en contenedores / poca RAM
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--window-size=1920,1080")

    # Silenciar logs de Chrome
    options.add_argument("--log-level=3")
    options.add_experimental_option("excludeSwitches", ["enable-logging"])

    driver = webdriver.Chrome(options=options)
    driver.get(url)

    try:
        # Esperar a que cargue el body
        element = WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.TAG_NAME, "body"))
        )
        info = element.text

        filename = f"peptide_{peptide_id}.txt"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(info)

        print(f"Información guardada en {filename}")
    except Exception as e:
        print(f"Error: {e}")
    finally:
        driver.quit()




