import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import re
import subprocess
import os
import threading
import time

class SignalAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Анализатор QPSK сигналов")
        
        # Настройки путей
        self.build_dir = "build"
        self.data_dir = "data"
        self.processing_app = "./data_processing"
        
        # Размеры
        self.main_window_size = "1600x1000"
        self.figure_sizes = [(10, 3), (14, 3)]
        self.frame_padding = {"padx": 10, "pady": 8}
        
        # Параметры демо-режима (QPSK)
        self.demo_params = {
            "fd": 1100.0,     # Sample freq
            "n": 10,        # Num info bits
            "vel": 10.0,    # Info velocity (bits/sec)
            "snr": 10.0,    # SNR
            "poly1": 37,    # Gold seq poly 1
            "poly2": 61     # Gold seq poly 2
        }
        
        # Параметры исследования BER
        self.research_params = {
            "fd": 1100.0,
            "n": 10,
            "vel": 10.0,
            "snr_min": 0.0,
            "snr_max": 20.0,
            "n_points": 10,
            "n_runs": 100,
            "poly1": 37,
            "poly2": 61
        }
        
        self.root.geometry(self.main_window_size)
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # Данные
        self.signal_iq = None
        self.correlations = None  # dict or list of 4 arrays
        self.in_bits = None
        self.out_bits = None
        self.ber_data = None
        
        self.current_mode = "demo"
        self.create_widgets()

    def create_widgets(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.demo_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.demo_frame, text="Демонстрация QPSK")
        
        self.research_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.research_frame, text="Исследование BER (QPSK)")
        
        self.create_demo_widgets()
        self.create_research_widgets()
        self.notebook.bind("<<NotebookTabChanged>>", self.on_tab_changed)

    def create_demo_widgets(self):
        # Параметры
        demo_params_frame = tk.LabelFrame(self.demo_frame, text="Параметры QPSK", font=('Arial', 10))
        demo_params_frame.pack(fill=tk.X, padx=10, pady=5)
        
        demo_param_labels = {
            "fd": "Частота дискретизации (fd):",
            "n": "Количество бит (n):",
            "vel": "Скорость (vel):",
            "snr": "SNR (дБ):",
            "poly1": "Полином 1 (poly1):",
            "poly2": "Полином 2 (poly2):"
        }
        
        self.demo_param_entries = {}
        row, col = 0, 0
        for key, label in demo_param_labels.items():
            frame = tk.Frame(demo_params_frame)
            frame.grid(row=row, column=col, padx=5, pady=2, sticky="w")
            tk.Label(frame, text=label, font=('Arial', 9)).pack(side=tk.LEFT)
            entry = tk.Entry(frame, width=10, font=('Arial', 9))
            entry.insert(0, str(self.demo_params[key]))
            entry.pack(side=tk.LEFT, padx=5)
            self.demo_param_entries[key] = entry
            
            col += 1
            if col > 2:
                col = 0
                row += 1
        
        # Кнопки
        btn_frame = tk.Frame(self.demo_frame)
        btn_frame.pack(pady=8)
        self.demo_start_btn = tk.Button(btn_frame, text="Старт", command=self.start_demo_processing,
                                       bg='green', fg='white', font=('Arial', 12))
        self.demo_start_btn.pack(side=tk.LEFT, padx=10)
        self.demo_exit_btn = tk.Button(btn_frame, text="Выход", command=self.on_closing,
                                      bg='red', fg='white', font=('Arial', 12))
        self.demo_exit_btn.pack(side=tk.LEFT, padx=10)
        
        # Биты: исходные и демодулированные
        bits_frame = tk.Frame(self.demo_frame)
        bits_frame.pack(fill=tk.X, padx=10, pady=5)
        
        tk.Label(bits_frame, text="Исходные биты:", font=('Arial', 10, 'bold')).pack(anchor='w')
        self.in_bits_text = tk.Text(bits_frame, height=2, width=100, font=('Courier', 10))
        self.in_bits_text.pack(fill=tk.X, pady=2)
        
        tk.Label(bits_frame, text="Демодулированные биты:", font=('Arial', 10, 'bold')).pack(anchor='w')
        self.out_bits_text = tk.Text(bits_frame, height=2, width=100, font=('Courier', 10))
        self.out_bits_text.pack(fill=tk.X, pady=2)
        
        # Графики
        plot_frame = tk.Frame(self.demo_frame)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        
        self.demo_frame1 = tk.Frame(plot_frame)
        self.demo_frame1.pack(fill=tk.BOTH, expand=True, **self.frame_padding)
        self.demo_frame2 = tk.Frame(plot_frame)
        self.demo_frame2.pack(fill=tk.BOTH, expand=True, **self.frame_padding)
        
        self.demo_fig1, self.demo_ax1 = plt.subplots(figsize=self.figure_sizes[0])
        self.demo_fig2, self.demo_ax2 = plt.subplots(figsize=self.figure_sizes[1])
        self.demo_fig1.tight_layout(pad=1.0)
        self.demo_fig2.tight_layout(pad=1.0)
        
        self.demo_canvas1 = FigureCanvasTkAgg(self.demo_fig1, master=self.demo_frame1)
        self.demo_canvas1.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.demo_canvas2 = FigureCanvasTkAgg(self.demo_fig2, master=self.demo_frame2)
        self.demo_canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self._create_status_bar()

    def create_research_widgets(self):
        research_params_frame = tk.LabelFrame(self.research_frame, text="Параметры BER исследования", font=('Arial', 10))
        research_params_frame.pack(fill=tk.X, padx=10, pady=5)
        
        research_param_labels = {
            "fd": "Частота дискретизации (fd):",
            "n": "Количество бит (n):",
            "vel": "Скорость (vel):",
            "snr_min": "Мин. SNR:",
            "snr_max": "Макс. SNR:",
            "n_points": "Точек:",
            "n_runs": "Испытаний:",
            "poly1": "Полином 1:",
            "poly2": "Полином 2:"
        }
        
        self.research_param_entries = {}
        row, col = 0, 0
        for key, label in research_param_labels.items():
            frame = tk.Frame(research_params_frame)
            frame.grid(row=row, column=col, padx=5, pady=2, sticky="w")
            tk.Label(frame, text=label, font=('Arial', 9)).pack(side=tk.LEFT)
            entry = tk.Entry(frame, width=10, font=('Arial', 9))
            entry.insert(0, str(self.research_params[key]))
            entry.pack(side=tk.LEFT, padx=5)
            self.research_param_entries[key] = entry
            
            col += 1
            if col > 2:
                col = 0
                row += 1
        
        btn_frame = tk.Frame(self.research_frame)
        btn_frame.pack(pady=8)
        self.research_start_btn = tk.Button(btn_frame, text="Начать исследование", 
                                           command=self.start_research_processing,
                                           bg='blue', fg='white', font=('Arial', 12))
        self.research_start_btn.pack(side=tk.LEFT, padx=10)
        self.research_exit_btn = tk.Button(btn_frame, text="Выход", command=self.on_closing,
                                          bg='red', fg='white', font=('Arial', 12))
        self.research_exit_btn.pack(side=tk.LEFT, padx=10)
        
        self.research_fig, self.research_ax = plt.subplots(figsize=(14, 8))
        self.research_fig.tight_layout(pad=2.0)
        self.research_canvas = FigureCanvasTkAgg(self.research_fig, master=self.research_frame)
        self.research_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self._create_status_bar()

    def _create_status_bar(self):
        if not hasattr(self, 'status_var'):
            self.status_var = tk.StringVar(value="Готов")
            self.status_bar = tk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
            self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
            self.progress_var = tk.DoubleVar()
            self.progress_bar = ttk.Progressbar(self.root, variable=self.progress_var, maximum=100)
            self.progress_bar.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=3)

    def on_tab_changed(self, event):
        self.current_mode = "demo" if self.notebook.index(self.notebook.select()) == 0 else "research"
        self.status_var.set(f"Режим: {'Демо QPSK' if self.current_mode == 'demo' else 'Исследование BER'}")

    def get_demo_parameters(self):
        try:
            return {
                "fd": float(self.demo_param_entries["fd"].get()),
                "n": int(self.demo_param_entries["n"].get()),
                "vel": float(self.demo_param_entries["vel"].get()),
                "snr": float(self.demo_param_entries["snr"].get()),
                "poly1": int(self.demo_param_entries["poly1"].get()),
                "poly2": int(self.demo_param_entries["poly2"].get())
            }
        except Exception as e:
            messagebox.showerror("Ошибка", f"Неверный параметр: {e}")
            return None

    def get_research_parameters(self):
        try:
            return {
                "fd": float(self.research_param_entries["fd"].get()),
                "n": int(self.research_param_entries["n"].get()),
                "vel": float(self.research_param_entries["vel"].get()),
                "snr_min": float(self.research_param_entries["snr_min"].get()),
                "snr_max": float(self.research_param_entries["snr_max"].get()),
                "n_points": int(self.research_param_entries["n_points"].get()),
                "n_runs": int(self.research_param_entries["n_runs"].get()),
                "poly1": int(self.research_param_entries["poly1"].get()),
                "poly2": int(self.research_param_entries["poly2"].get())
            }
        except Exception as e:
            messagebox.showerror("Ошибка", f"Неверный параметр: {e}")
            return None

    def check_directories(self):
        if not os.path.exists(self.build_dir):
            messagebox.showerror("Ошибка", f"Нет директории {self.build_dir}")
            return False
        if not os.path.exists(os.path.join(self.build_dir, self.processing_app)):
            messagebox.showerror("Ошибка", f"Нет {self.processing_app} в {self.build_dir}")
            return False
        os.makedirs(self.data_dir, exist_ok=True)
        return True

    def run_data_processing(self, args):
        try:
            orig = os.getcwd()
            os.chdir(self.build_dir)
            result = subprocess.run(args, capture_output=True, text=True, timeout=60)
            os.chdir(orig)
            if result.returncode != 0:
                messagebox.showerror("Ошибка", result.stderr)
                return False
            return True
        except Exception as e:
            messagebox.showerror("Ошибка", str(e))
            if 'orig' in locals():
                os.chdir(orig)
            return False

    def run_demo_processing(self, params):
        args = [
            self.processing_app,
            str(params["fd"]),
            str(params["n"]),
            str(params["vel"]),
            str(params["snr"]),
            str(params["poly1"]),
            str(params["poly2"])
        ]
        return self.run_data_processing(args)

    def run_research_point(self, params, snr, idx, total):
        args = [
            self.processing_app,
            str(params["fd"]),
            str(params["n"]),
            str(params["vel"]),
            str(snr),
            str(params["poly1"]),
            str(params["poly2"]),
            str(params["n_runs"])
        ]
        self.status_var.set(f"Точка {idx+1}/{total}, SNR={snr:.2f} дБ")
        return self.run_data_processing(args)

    def parse_real_flat_file(self, filename):
        """Читает файл с чередующимися I,Q,I,Q... как вещественные числа → complex array"""
        try:
            data = np.loadtxt(filename, delimiter=',', dtype=float)
            if len(data) % 2 != 0:
                data = data[:-1]  # drop last if odd
            i_vals = data[0::2]
            q_vals = data[1::2]
            return i_vals + 1j * q_vals
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка чтения {filename}: {e}")
            return None

    def parse_real_txt_file(self, filename):
        """Обычный парсинг вещественных чисел (для корреляций и BER)"""
        try:
            data = np.loadtxt(filename, delimiter=',', dtype=float)
            return data if data.ndim == 1 else data.flatten()
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка чтения {filename}: {e}")
            return None

    def parse_bits_file(self, filename):
        """Читает биты как целые 0/1"""
        try:
            with open(filename, 'r') as f:
                content = f.read().replace(',', ' ').replace('(', '').replace(')', '')
                bits = [int(x) for x in content.split() if x.strip() in ('0', '1')]
            return bits
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка чтения битов {filename}: {e}")
            return None

    def load_demo_files(self):
        files = {
            "signalIQ.txt": "signal",
            "correlation1.txt": "corr1",
            "correlation2.txt": "corr2",
            "correlation3.txt": "corr3",
            "correlation4.txt": "corr4",
            "in_bits.txt": "in_bits",
            "out_bits.txt": "out_bits"
        }
        
        for fname, key in files.items():
            path = os.path.join(self.data_dir, fname)
            if not os.path.exists(path):
                messagebox.showerror("Ошибка", f"Файл отсутствует: {fname}")
                return False
            
            if key == "signal":
                self.signal_iq = self.parse_real_flat_file(path)
                if self.signal_iq is None: return False
            elif key in ["corr1", "corr2", "corr3", "corr4"]:
                corr = self.parse_real_txt_file(path)
                if corr is None: return False
                if self.correlations is None:
                    self.correlations = []
                self.correlations.append(corr)
            elif key == "in_bits":
                self.in_bits = self.parse_bits_file(path)
                if self.in_bits is None: return False
            elif key == "out_bits":
                self.out_bits = self.parse_bits_file(path)
                if self.out_bits is None: return False
        
        return True

    def load_ber_file(self):
        path = os.path.join(self.data_dir, "ber_qpsk.txt")
        if not os.path.exists(path):
            messagebox.showerror("Ошибка", "Файл ber_qpsk.txt не найден")
            return None
        return self.parse_real_txt_file(path)

    def show_demo_plots(self):
        # Очистка
        self.demo_ax1.clear()
        self.demo_ax2.clear()
        
        # 1. I/Q сигнал
        time = np.arange(len(self.signal_iq))
        self.demo_ax1.plot(time, self.signal_iq.real, 'b-', label='I (In-phase)', linewidth=1)
        self.demo_ax1.plot(time, self.signal_iq.imag, 'r-', label='Q (Quadrature)', linewidth=1)
        self.demo_ax1.set_title('Зашумленный QPSK сигнал (I/Q компоненты)', fontsize=12)
        self.demo_ax1.set_xlabel('Отсчеты')
        self.demo_ax1.set_ylabel('Амплитуда')
        self.demo_ax1.legend()
        self.demo_ax1.grid(True, alpha=0.3)
        
        # 2. 4 корреляции
        labels = ["Символ '00'", "Символ '01'", "Символ '10'", "Символ '11'"]
        colors = ['blue', 'red', 'green', 'purple']
        for i, (corr, label, color) in enumerate(zip(self.correlations, labels, colors)):
            self.demo_ax2.plot(corr, color=color, label=label, linewidth=1.5)
        self.demo_ax2.set_title('Отклики 4-х корреляторов', fontsize=12)
        self.demo_ax2.set_xlabel('Отсчеты')
        self.demo_ax2.set_ylabel('Амплитуда корреляции')
        self.demo_ax2.legend()
        self.demo_ax2.grid(True, alpha=0.3)
        
        self.demo_fig1.tight_layout()
        self.demo_fig2.tight_layout()
        self.demo_canvas1.draw()
        self.demo_canvas2.draw()
        
        # Биты
        self.in_bits_text.delete(1.0, tk.END)
        self.in_bits_text.insert(tk.END, ''.join(map(str, self.in_bits)) if self.in_bits else "")
        self.out_bits_text.delete(1.0, tk.END)
        self.out_bits_text.insert(tk.END, ''.join(map(str, self.out_bits)) if self.out_bits else "")

    def show_research_plot(self, ber_data, snr_values):
        self.research_ax.clear()
        self.research_ax.semilogy(snr_values, ber_data, 'bo-', linewidth=2, markersize=6)
        self.research_ax.set_title('SER (Symbol Error Ratio) для QPSK модуляции', fontsize=14)
        self.research_ax.set_xlabel('SNR (дБ)')
        self.research_ax.set_ylabel('Symbol Error Ratio (SER)')
        self.research_ax.grid(True, which='both', alpha=0.5)
        self.research_ax.set_yscale('log')
        self.research_fig.tight_layout()
        self.research_canvas.draw()

    def demo_processing_thread(self):
        self.signal_iq = None
        self.correlations = []  # Используйте список, а не None
        self.in_bits = None
        self.out_bits = None

        # ОЧИСТКА ГРАФИКОВ
        self.demo_ax1.clear()
        self.demo_ax2.clear()
        self.demo_canvas1.draw()
        self.demo_canvas2.draw()
        try:
            self.demo_start_btn.config(state=tk.DISABLED)
            params = self.get_demo_parameters()
            if not params: return
            if not self.check_directories(): return
            
            if not self.run_demo_processing(params): return
            if not self.load_demo_files(): return
            self.show_demo_plots()
            messagebox.showinfo("Успех", "Демонстрация завершена!")
        except Exception as e:
            messagebox.showerror("Ошибка", str(e))
        finally:
            self.demo_start_btn.config(state=tk.NORMAL)

    def research_processing_thread(self):
        try:
            self.research_start_btn.config(state=tk.DISABLED)
            params = self.get_research_parameters()
            if not params: return
            if params["n_points"] <= 0 or params["snr_min"] >= params["snr_max"]:
                messagebox.showerror("Ошибка", "Неверные параметры SNR или точек")
                return
            if not self.check_directories(): return
            
            # Удалим старый файл
            ber_path = os.path.join(self.data_dir, "ber_qpsk.txt")
            if os.path.exists(ber_path):
                os.remove(ber_path)
            
            snr_vals = np.linspace(params["snr_min"], params["snr_max"], params["n_points"])
            total = len(snr_vals)
            
            for i, snr in enumerate(snr_vals):
                self.progress_var.set(5 + 90 * i / total)
                if not self.run_research_point(params, snr, i, total):
                    return
                time.sleep(0.05)
            
            ber = self.load_ber_file()
            if ber is None: return
            if len(ber) != len(snr_vals):
                messagebox.showerror("Ошибка", "Несовпадение длины BER и SNR")
                return
            
            self.show_research_plot(ber, snr_vals)
            messagebox.showinfo("Успех", "Исследование BER завершено!")
        except Exception as e:
            messagebox.showerror("Ошибка", str(e))
        finally:
            self.research_start_btn.config(state=tk.NORMAL)

    def start_demo_processing(self):
        threading.Thread(target=self.demo_processing_thread, daemon=True).start()

    def start_research_processing(self):
        threading.Thread(target=self.research_processing_thread, daemon=True).start()

    def on_closing(self):
        if messagebox.askokcancel("Выход", "Завершить работу?"):
            self.root.quit()  # This stops the mainloop
            self.root.destroy()  # This destroys the window
            # Exit the Python process
            import sys
            sys.exit(0)

def main():
    root = tk.Tk()
    app = SignalAnalyzerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()  