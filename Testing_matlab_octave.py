# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:27:57 2024

@author: diego
"""

import subprocess

# Define el comando para ejecutar Octave con tu script de MATLAB
comando_octave = 'octave --eval "mi_script_matlab.m"'

# Ejecuta el comando en la línea de comandos y captura la salida
resultado = subprocess.check_output(comando_octave, shell=True)

# Convierte la salida a una cadena legible en Python si es necesario
resultado_decodificado = resultado.decode('utf-8')

# Imprime o maneja la salida de Octave en Python
print(resultado_decodificado)