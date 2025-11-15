import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import re
import cmath
import subprocess
import os
import threading
import time

class SignalAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Анализатор сигналов - Демо и Исследование")
        
        # Настройки путей
        self.build_dir = "build"
        self.data_dir = "data"
        self.processing_app = "./data_processing"
        
        # Настройки размеров
        self.main_window_size = "1600x1000"
        self.figure_sizes = [(14, 2), (14, 2), (14, 2)]
        self.frame_padding = {"padx": 13, "pady": 8}
        
        # Параметры по умолчанию для демо режима
        self.demo_params = {
            "fd": 20.0,    # Sample freq
            "f": 10.0,     # Carrier freq
            "n": 10,       # Num info bits
            "vel": 10.0,   # Info velocity
            "dt": 0.0,     # Time offset
            "snr1": 10.0,  # SNR for signal 1
            "snr2": 10.0,  # SNR for signal 2
            "type": 0      # modulation type 
        }
        
        # Параметры по умолчанию для режима исследования
        self.research_params = {
            "fd": 20.0,        # Sample freq
            "f": 10.0,         # Carrier freq
            "n": 100,          # Num info bits
            "vel": 10.0,       # Info velocity
            "snr_static": 10.0,# Статическое SNR для референсного сигнала
            "snr_min": 0.0,    # Минимальное SNR для BER
            "snr_max": 20.0,   # Максимальное SNR для BER
            "n_points": 10,    # Количество точек на графике
            "n_runs": 100      # Количество испытаний на точку
        }
        
        self.root.geometry(self.main_window_size)
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        self.signal1_data = None
        self.signal2_data = None
        self.correlation_data = None
        self.ber_data = None
        
        self.current_mode = "demo"  # "demo" или "research"
        
        self.create_widgets()
        
    def create_widgets(self):
        # Создаем Notebook для вкладок
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Вкладка для демо режима
        self.demo_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.demo_frame, text="Демонстрационный режим")
        
        # Вкладка для режима исследования
        self.research_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.research_frame, text="Режим исследования")
        
        # Создаем виджеты для демо режима
        self.create_demo_widgets()
        
        # Создаем виджеты для режима исследования
        self.create_research_widgets()
        
        # Биндим событие переключения вкладок
        self.notebook.bind("<<NotebookTabChanged>>", self.on_tab_changed)
    
    def create_demo_widgets(self):
        # Фрейм для параметров демо режима
        demo_params_frame = tk.LabelFrame(self.demo_frame, text="Параметры обработки", font=('Arial', 10))
        demo_params_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Создаем поля для параметров демо режима
        demo_param_labels = {
            "fd": "Частота дискретизации (fd):",
            "f": "Несущая частота (f):",
            "n": "Количество бит информации (n):",
            "vel": "Скорость информации (vel):",
            "dt": "Временное смещение (dt):",
            "snr1": "SNR для сигнала 1 (snr1):",
            "snr2": "SNR для сигнала 2 (snr2):",
            "type": "Тип модуляции (0 - АМ, 1 - ФМ-2, 2 - МЧМ):"
        }
        
        self.demo_param_entries = {}
        row = 0
        col = 0
        
        for i, (key, label) in enumerate(demo_param_labels.items()):
            frame = tk.Frame(demo_params_frame)
            frame.grid(row=row, column=col, padx=5, pady=2, sticky="w")
            
            tk.Label(frame, text=label, font=('Arial', 9)).pack(side=tk.LEFT)
            entry = tk.Entry(frame, width=10, font=('Arial', 9))
            entry.insert(0, str(self.demo_params[key]))
            entry.pack(side=tk.LEFT, padx=5)
            self.demo_param_entries[key] = entry
            
            col += 1
            if col > 2:  # 3 колонки в ряду
                col = 0
                row += 1
        
        # Фрейм для кнопок демо режима
        demo_button_frame = tk.Frame(self.demo_frame)
        demo_button_frame.pack(pady=10)
        
        # Кнопка старт для демо режима
        self.demo_start_btn = tk.Button(demo_button_frame, text="Старт обработки", 
                                       command=self.start_demo_processing,
                                       bg='green', fg='white', font=('Arial', 12))
        self.demo_start_btn.pack(side=tk.LEFT, padx=10)
        
        # Кнопка выхода
        self.demo_exit_btn = tk.Button(demo_button_frame, text="Выход", 
                                      command=self.on_closing,
                                      bg='red', fg='white', font=('Arial', 12))
        self.demo_exit_btn.pack(side=tk.LEFT, padx=10)
        
        # Создаем фреймы для трех отдельных графиков в демо режиме
        self.demo_plot_frame = tk.Frame(self.demo_frame)
        self.demo_plot_frame.pack(fill=tk.BOTH, expand=True)
        
        # Фреймы для каждого графика
        self.demo_frame1 = tk.Frame(self.demo_plot_frame)
        self.demo_frame1.pack(side=tk.TOP, fill=tk.BOTH, expand=True, **self.frame_padding)
        
        self.demo_frame2 = tk.Frame(self.demo_plot_frame)
        self.demo_frame2.pack(side=tk.TOP, fill=tk.BOTH, expand=True, **self.frame_padding)
        
        self.demo_frame3 = tk.Frame(self.demo_plot_frame)
        self.demo_frame3.pack(side=tk.TOP, fill=tk.BOTH, expand=True, **self.frame_padding)
        
        # Создаем отдельные фигуры для демо режима
        self.demo_fig1, self.demo_ax1 = plt.subplots(figsize=self.figure_sizes[0])
        self.demo_fig2, self.demo_ax2 = plt.subplots(figsize=self.figure_sizes[1])
        self.demo_fig3, self.demo_ax3 = plt.subplots(figsize=self.figure_sizes[2])
        
        # Настраиваем отступы внутри фигур
        figure_padding = 1.0
        self.demo_fig1.tight_layout(pad=figure_padding)
        self.demo_fig2.tight_layout(pad=figure_padding)
        self.demo_fig3.tight_layout(pad=figure_padding)
        
        # Canvas для каждого графика в демо режиме
        self.demo_canvas1 = FigureCanvasTkAgg(self.demo_fig1, master=self.demo_frame1)
        self.demo_canvas1.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.demo_canvas2 = FigureCanvasTkAgg(self.demo_fig2, master=self.demo_frame2)
        self.demo_canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.demo_canvas3 = FigureCanvasTkAgg(self.demo_fig3, master=self.demo_frame3)
        self.demo_canvas3.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def create_research_widgets(self):
        # Фрейм для параметров исследования
        research_params_frame = tk.LabelFrame(self.research_frame, text="Параметры исследования BER", font=('Arial', 10))
        research_params_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Создаем поля для параметров исследования
        research_param_labels = {
            "fd": "Частота дискретизации (fd):",
            "f": "Несущая частота (f):",
            "n": "Количество бит информации (n):",
            "vel": "Скорость информации (vel):",
            "snr_static": "Статическое SNR (snr_static):",
            "snr_min": "Минимальное SNR (snr_min):",
            "snr_max": "Максимальное SNR (snr_max):",
            "n_points": "Количество точек (n_points):",
            "n_runs": "Количество испытаний (n_runs):"
        }
        
        self.research_param_entries = {}
        row = 0
        col = 0
        
        for i, (key, label) in enumerate(research_param_labels.items()):
            frame = tk.Frame(research_params_frame)
            frame.grid(row=row, column=col, padx=5, pady=2, sticky="w")
            
            tk.Label(frame, text=label, font=('Arial', 9)).pack(side=tk.LEFT)
            entry = tk.Entry(frame, width=10, font=('Arial', 9))
            entry.insert(0, str(self.research_params[key]))
            entry.pack(side=tk.LEFT, padx=5)
            self.research_param_entries[key] = entry
            
            col += 1
            if col > 2:  # 3 колонки в ряду
                col = 0
                row += 1
        
        # Фрейм для кнопок исследования
        research_button_frame = tk.Frame(self.research_frame)
        research_button_frame.pack(pady=10)
        
        # Кнопка начала исследования
        self.research_start_btn = tk.Button(research_button_frame, text="Начать исследование", 
                                           command=self.start_research_processing,
                                           bg='blue', fg='white', font=('Arial', 12))
        self.research_start_btn.pack(side=tk.LEFT, padx=10)
        
        # Кнопка выхода
        self.research_exit_btn = tk.Button(research_button_frame, text="Выход", 
                                          command=self.on_closing,
                                          bg='red', fg='white', font=('Arial', 12))
        self.research_exit_btn.pack(side=tk.LEFT, padx=10)
        
        # Создаем фрейм для графика BER в режиме исследования
        self.research_plot_frame = tk.Frame(self.research_frame)
        self.research_plot_frame.pack(fill=tk.BOTH, expand=True)
        
        # Создаем фигуру для графика BER
        self.research_fig, self.research_ax = plt.subplots(figsize=(14, 8))
        self.research_fig.tight_layout(pad=3.0)
        
        # Canvas для графика BER
        self.research_canvas = FigureCanvasTkAgg(self.research_fig, master=self.research_plot_frame)
        self.research_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Статус бар и прогресс бар (общие для обоих режимов)
        self.status_var = tk.StringVar()
        self.status_var.set("Готов к работе")
        self.status_bar = tk.Label(self.root, textvariable=self.status_var, 
                                  relief=tk.SUNKEN, anchor=tk.W, font=('Arial', 10))
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(self.root, variable=self.progress_var, 
                                          maximum=100, mode='determinate')
        self.progress_bar.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
    
    def on_tab_changed(self, event):
        """Обработчик переключения вкладок"""
        selected_tab = self.notebook.index(self.notebook.select())
        if selected_tab == 0:
            self.current_mode = "demo"
            self.status_var.set("Режим: Демонстрационный")
        else:
            self.current_mode = "research"
            self.status_var.set("Режим: Исследование BER")
    
    def get_demo_parameters(self):
        """Получение параметров демо режима из полей ввода"""
        try:
            params = {}
            for key, entry in self.demo_param_entries.items():
                value = entry.get().strip()
                if key == "n" or key == "type":
                    params[key] = int(value) if value else 0
                else:
                    params[key] = float(value) if value else 0.0
            return params
        except ValueError as e:
            messagebox.showerror("Ошибка", f"Некорректное значение параметра: {str(e)}")
            return None
    
    def get_research_parameters(self):
        """Получение параметров исследования из полей ввода"""
        try:
            params = {}
            for key, entry in self.research_param_entries.items():
                value = entry.get().strip()
                if key in ["n", "n_points", "n_runs"]:
                    params[key] = int(value) if value else 0
                else:
                    params[key] = float(value) if value else 0.0
            return params
        except ValueError as e:
            messagebox.showerror("Ошибка", f"Некорректное значение параметра: {str(e)}")
            return None
    
    def on_closing(self):
        """Обработчик закрытия приложения"""
        if messagebox.askokcancel("Выход", "Вы уверены, что хотите выйти?"):
            self.root.quit()
            self.root.destroy()
    
    def check_directories(self):
        """Проверка существования необходимых директорий"""
        if not os.path.exists(self.build_dir):
            messagebox.showerror("Ошибка", f"Директория {self.build_dir} не существует!")
            return False
        
        if not os.path.exists(os.path.join(self.build_dir, self.processing_app)):
            messagebox.showerror("Ошибка", 
                               f"Утилита {self.processing_app} не найдена в {self.build_dir}!")
            return False
        
        # Создаем data директорию если ее нет
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
            self.status_var.set(f"Создана директория {self.data_dir}")
        
        return True
    
    def cleanup_ber_files(self):
        """Очистка старых BER файлов перед началом исследования"""
        ber_files = ["ber_am.txt", "ber_fm.txt", "ber_pm.txt"]
        
        for filename in ber_files:
            filepath = os.path.join(self.data_dir, filename)
            if os.path.exists(filepath):
                try:
                    os.remove(filepath)
                    self.status_var.set(f"Удален старый файл: {filename}")
                except Exception as e:
                    messagebox.showwarning("Предупреждение", 
                                         f"Не удалось удалить файл {filename}: {str(e)}")
    
    def run_data_processing(self, args):
        """Запуск внешней утилиты обработки данных с аргументами"""
        try:
            # Переходим в build директорию
            original_dir = os.getcwd()
            os.chdir(self.build_dir)
            
            # Запускаем утилиту
            result = subprocess.run(args, 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=30)  # Таймаут 30 секунд
            
            # Возвращаемся обратно
            os.chdir(original_dir)
            
            if result.returncode != 0:
                error_msg = f"Ошибка выполнения {self.processing_app}:\n{result.stderr}"
                messagebox.showerror("Ошибка", error_msg)
                return False
            
            return True
            
        except subprocess.TimeoutExpired:
            os.chdir(original_dir)
            messagebox.showerror("Ошибка", "Утилита превысила время выполнения!")
            return False
        except Exception as e:
            os.chdir(original_dir)
            messagebox.showerror("Ошибка", f"Ошибка при запуске утилиты: {str(e)}")
            return False
    
    def run_demo_processing(self, params):
        """Запуск обработки для демо режима"""
        # Формируем аргументы командной строки
        args = [
            self.processing_app,
            str(params["fd"]),
            str(params["f"]),
            str(params["n"]),
            str(params["vel"]),
            str(params["dt"]),
            str(params["snr1"]),
            str(params["snr2"]),
            str(params["type"])
        ]
        
        self.status_var.set("Запуск data_processing...")
        self.progress_var.set(20)
        
        return self.run_data_processing(args)
    
    def run_research_point(self, params, snr2, point_index, total_points):
        """Запуск одной точки исследования для всех типов модуляции"""
        # Формируем аргументы командной строки для исследования
        # dt и type не указываются - программа сама генерирует данные для всех модуляций
        args = [
            self.processing_app,
            str(params["fd"]),
            str(params["f"]),
            str(params["n"]),
            str(params["vel"]),
            str(params["snr_static"]),  # snr1 - статическое значение
            str(snr2),                  # snr2 - переменное значение
            str(params["n_runs"])       # Количество испытаний
        ]
        
        self.status_var.set(f"Точка {point_index+1}/{total_points}, SNR={snr2:.2f} дБ")
        
        return self.run_data_processing(args)
    
    def parse_complex_txt_file(self, filename):
        """Парсинг файла с комплексными числами в формате (real,imag)"""
        complex_data = []
        try:
            with open(filename, 'r', encoding='utf-8') as file:
                content = file.read()
                pattern = r'\(([-\d\.]+),([-\d\.]+)\)'
                matches = re.findall(pattern, content)
                
                for real_str, imag_str in matches:
                    try:
                        real = float(real_str.strip())
                        imag = float(imag_str.strip())
                        complex_data.append(complex(real, imag))
                    except ValueError:
                        continue
            
            return np.array(complex_data)
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при чтении файла {filename}:\n{str(e)}")
            return None
    
    def parse_real_txt_file(self, filename):
        """Парсинг файла с действительными числами"""
        real_data = []
        try:
            with open(filename, 'r', encoding='utf-8') as file:
                content = file.read()
                numbers = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', content)
                
                for num_str in numbers:
                    try:
                        real_data.append(float(num_str))
                    except ValueError:
                        continue
            
            return np.array(real_data)
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при чтении файла {filename}:\n{str(e)}")
            return None
    
    def load_demo_files(self):
        """Загрузка сгенерированных файлов для демо режима"""
        files_to_load = [
            ("first_data.txt", "signal1"),
            ("second_data.txt", "signal2"), 
            ("correlation.txt", "correlation")
        ]
        
        for filename, data_type in files_to_load:
            filepath = os.path.join(self.data_dir, filename)
            
            if not os.path.exists(filepath):
                messagebox.showerror("Ошибка", f"Файл {filename} не найден в {self.data_dir}!")
                return False
            
            self.status_var.set(f"Загрузка {filename}...")
            
            if data_type in ["signal1", "signal2"]:
                data = self.parse_complex_txt_file(filepath)
            else:
                data = self.parse_real_txt_file(filepath)
            
            if data is None or len(data) == 0:
                messagebox.showerror("Ошибка", f"Не удалось загрузить данные из {filename}!")
                return False
            
            if data_type == "signal1":
                self.signal1_data = data
            elif data_type == "signal2":
                self.signal2_data = data
            else:
                self.correlation_data = data
            
            self.progress_var.set(self.progress_var.get() + 15)
            time.sleep(0.1)
        
        return True
    
    def load_ber_files(self):
        """Загрузка BER данных для исследования"""
        ber_files = ["ber_am.txt", "ber_fm.txt", "ber_pm.txt"]
        ber_data = {}
        
        for filename in ber_files:
            filepath = os.path.join(self.data_dir, filename)
            
            if not os.path.exists(filepath):
                messagebox.showerror("Ошибка", f"Файл {filename} не найден в {self.data_dir}!")
                return None
            
            data = self.parse_real_txt_file(filepath)
            if data is None or len(data) == 0:
                messagebox.showerror("Ошибка", f"Не удалось загрузить данные из {filename}!")
                return None
            
            # Извлекаем название модуляции из имени файла
            mod_name = filename.replace("ber_", "").replace(".txt", "").upper()
            ber_data[mod_name] = data
        
        return ber_data
    
    def show_demo_plots(self):
        """Отображение графиков для демо режима"""
        try:
            # Очищаем графики
            self.demo_ax1.clear()
            self.demo_ax2.clear()
            self.demo_ax3.clear()
            
            # 1. График первого комплексного сигнала
            time1 = np.arange(len(self.signal1_data))
            self.demo_ax1.plot(time1, np.real(self.signal1_data), 'b-', label='Сигнал 1 (I)', linewidth=1)
            self.demo_ax1.plot(time1, np.imag(self.signal1_data), 'r-', label='Сигнал 1 (Q)', linewidth=1)
            self.demo_ax1.set_title('Сигнал 1 - Комплексные I/Q компоненты', fontsize=12)
            self.demo_ax1.set_xlabel('Отсчеты', fontsize=10)
            self.demo_ax1.set_ylabel('Амплитуда', fontsize=10)
            self.demo_ax1.legend(fontsize=9)
            self.demo_ax1.grid(True, alpha=0.3)
            
            # 2. График второго комплексного сигнала
            time2 = np.arange(len(self.signal2_data))
            self.demo_ax2.plot(time2, np.real(self.signal2_data), 'g-', label='Сигнал 2 (I)', linewidth=1)
            self.demo_ax2.plot(time2, np.imag(self.signal2_data), 'm-', label='Сигнал 2 (Q)', linewidth=1)
            self.demo_ax2.set_title('Сигнал 2 - Комплексные I/Q компоненты', fontsize=12)
            self.demo_ax2.set_xlabel('Отсчеты', fontsize=10)
            self.demo_ax2.set_ylabel('Амплитуда', fontsize=10)
            self.demo_ax2.legend(fontsize=9)
            self.demo_ax2.grid(True, alpha=0.3)
            
            # 3. График модуля корреляции
            time_corr = np.arange(len(self.correlation_data))
            self.demo_ax3.plot(time_corr, self.correlation_data, 'orange', linewidth=2)
            self.demo_ax3.set_title('Модуль корреляции сигналов', fontsize=12)
            self.demo_ax3.set_xlabel('Отсчеты', fontsize=10)
            self.demo_ax3.set_ylabel('Амплитуда корреляции', fontsize=10)
            self.demo_ax3.grid(True, alpha=0.3)
            
            # Помечаем максимум корреляции
            max_idx = np.argmax(self.correlation_data)
            max_val = self.correlation_data[max_idx]
            self.demo_ax3.plot(max_idx, max_val, 'ro', markersize=8)
            self.demo_ax3.annotate(f'Максимум: {max_val:.3f}\nОтсчет: {max_idx}', 
                                 xy=(max_idx, max_val),
                                 xytext=(max_idx + len(self.correlation_data)*0.1, max_val*0.9),
                                 arrowprops=dict(arrowstyle='->', color='red'),
                                 bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.3),
                                 fontsize=9)
            
            # Обновляем canvas
            self.demo_fig1.tight_layout()
            self.demo_fig2.tight_layout()
            self.demo_fig3.tight_layout()
            
            self.demo_canvas1.draw()
            self.demo_canvas2.draw()
            self.demo_canvas3.draw()
            
            self.status_var.set("Графики успешно отображены")
            self.progress_var.set(100)
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при отображении графиков:\n{str(e)}")
    
    def show_research_plot(self, ber_data, snr_values):
        """Отображение графика BER для исследования"""
        try:
            # Очищаем график
            self.research_ax.clear()
            
            # Цвета для разных типов модуляции
            colors = {'AM': 'blue', 'FM': 'red', 'PM': 'green'}
            markers = {'AM': 'o', 'FM': 's', 'PM': '^'}
            labels = {'AM': 'Амплитудная модуляция (АМ)', 
                     'FM': 'Частотная модуляция (МЧМ)', 
                     'PM': 'Фазовая модуляция (ФМ-2)'}
            
            # Строим кривые для каждого типа модуляции
            for mod_type, ber_values in ber_data.items():
                if len(ber_values) == len(snr_values):
                    self.research_ax.semilogy(snr_values, ber_values, 
                                            color=colors.get(mod_type, 'black'),
                                            marker=markers.get(mod_type, 'o'),
                                            label=labels.get(mod_type, mod_type),
                                            linewidth=2,
                                            markersize=6)
            
            self.research_ax.set_title('Зависимость BER от SNR для различных типов модуляции', fontsize=14)
            self.research_ax.set_xlabel('SNR (дБ)', fontsize=12)
            self.research_ax.set_ylabel('Bit Error Rate (BER)', fontsize=12)
            self.research_ax.legend(fontsize=10)
            self.research_ax.grid(True, alpha=0.3, which='both')
            self.research_ax.set_yscale('log')
            
            # Обновляем canvas
            self.research_fig.tight_layout()
            self.research_canvas.draw()
            
            self.status_var.set("График BER успешно отображен")
            self.progress_var.set(100)
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при отображении графика BER:\n{str(e)}")
    
    def demo_processing_thread(self):
        """Поток обработки данных для демо режима"""
        try:
            self.demo_start_btn.config(state=tk.DISABLED)
            self.progress_var.set(0)
            
            # Получаем параметры
            self.status_var.set("Чтение параметров...")
            params = self.get_demo_parameters()
            if params is None:
                return
            
            # Шаг 1: Проверка директорий
            self.status_var.set("Проверка директорий...")
            if not self.check_directories():
                return
            
            self.progress_var.set(10)
            
            # Шаг 2: Запуск внешней утилиты с параметрами
            self.status_var.set("Запуск data_processing с параметрами...")
            if not self.run_demo_processing(params):
                return
            
            # Шаг 3: Загрузка сгенерированных файлов
            self.status_var.set("Загрузка данных...")
            if not self.load_demo_files():
                return
            
            # Шаг 4: Отображение графиков
            self.status_var.set("Построение графиков...")
            self.show_demo_plots()
            
            messagebox.showinfo("Успех", "Обработка данных завершена успешно!")
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка в процессе обработки: {str(e)}")
        finally:
            self.demo_start_btn.config(state=tk.NORMAL)
    
    def research_processing_thread(self):
        """Поток обработки данных для исследования"""
        try:
            self.research_start_btn.config(state=tk.DISABLED)
            self.progress_var.set(0)
            
            # Получаем параметры
            self.status_var.set("Чтение параметров исследования...")
            params = self.get_research_parameters()
            if params is None:
                return
            
            # Проверяем параметры
            if params["n_points"] <= 0:
                messagebox.showerror("Ошибка", "Количество точек должно быть больше 0")
                return
            
            if params["snr_min"] >= params["snr_max"]:
                messagebox.showerror("Ошибка", "Минимальное SNR должно быть меньше максимального")
                return
            
            # Шаг 1: Проверка директорий
            self.status_var.set("Проверка директорий...")
            if not self.check_directories():
                return
            
            # Шаг 2: Очистка старых BER файлов
            self.status_var.set("Очистка старых данных BER...")
            self.cleanup_ber_files()
            
            self.progress_var.set(5)
            
            # Шаг 3: Генерируем значения SNR
            snr_step = (params["snr_max"] - params["snr_min"]) / params["n_points"]
            snr_values = [params["snr_min"] + i * snr_step for i in range(params["n_points"])]
            
            total_points = len(snr_values)
            
            # Шаг 4: Итеративный запуск утилиты для каждой точки SNR
            # Программа сама генерирует данные для всех трех типов модуляции
            for i, snr2 in enumerate(snr_values):
                point_progress = (i / total_points) * 90 + 5  # 5-95%
                self.progress_var.set(point_progress)
                
                if not self.run_research_point(params, snr2, i, total_points):
                    return
                
                time.sleep(0.1)  # Небольшая пауза между запусками
            
            # Шаг 5: Загрузка BER данных
            self.status_var.set("Загрузка данных BER...")
            self.progress_var.set(95)
            
            ber_data = self.load_ber_files()
            if ber_data is None:
                return
            
            # Шаг 6: Отображение графика
            self.status_var.set("Построение графика BER...")
            self.show_research_plot(ber_data, snr_values)
            
            messagebox.showinfo("Успех", "Исследование BER завершено успешно!")
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка в процессе исследования: {str(e)}")
        finally:
            self.research_start_btn.config(state=tk.NORMAL)
    
    def start_demo_processing(self):
        """Запуск процесса обработки в отдельном потоке для демо режима"""
        thread = threading.Thread(target=self.demo_processing_thread)
        thread.daemon = True
        thread.start()
    
    def start_research_processing(self):
        """Запуск процесса исследования в отдельном потоке"""
        thread = threading.Thread(target=self.research_processing_thread)
        thread.daemon = True
        thread.start()

def main():
    root = tk.Tk()
    app = SignalAnalyzerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()