# Modelo COVID-19
Este repositorio contiene un modelo matemático para el ajuste de datos epidemiológicos de COVID-19 en Colima-Villa de Álvarez.

📌 Descripción
El modelo utiliza métodos estadísticos avanzados para estimar la evolución de la epidemia, considerando distintos estados de la enfermedad y medidas de intervención.

📊 Suposiciones del modelo
Ajuste a datos de infecciosos, aislados (ambulatorios y hospitalizados), recuperados, defunciones y contagios acumulados.
Se asume que todas las personas con diagnóstico confirmado se aíslan.
Se consideran infecciosas a las personas desde el primer día de síntomas.
Se modelan explícitamente los casos ambulatorios y hospitalizados usando distribuciones de tiempo gamma.
Se estima que una proporción de los casos no son confirmados.
Se emplea el método de Cadenas de Markov y Monte Carlo (MCMC) para la estimación de parámetros.

🔢 Parámetros principales
Algunos de los parámetros clave del modelo incluyen:

Tasa de contagio (beta): 0.35 - 0.5
Factor de confirmación (factorConf): 64 - 70
Proporción de casos ambulatorios (pA): 0.8 - 0.86
Tasa de progresión a hospitalización (sigma): 0.26 - 0.27
Tasa de recuperación en aislados ambulatorios (gammaAC): 0.14 - 0.16
Tasa de recuperación en hospitalizados (gammaHC): 0.098 - 0.102
Mortalidad hospitalaria (muC): 0.2 - 0.35

⚙️ Metodología
El modelo se ajusta a los datos reales mediante la optimización de parámetros usando MCMC. Se tienen en cuenta intervenciones gubernamentales y cambios en la tasa de contagio.
