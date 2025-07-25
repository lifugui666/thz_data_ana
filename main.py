import sys
from PyQt5.QtWidgets import QMainWindow,QApplication,QWidget
from Ui_main_windows import Ui_MainWindow  #导入你写的界面类
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QVBoxLayout, QLabel
from PyQt5.QtGui import QImage, QPixmap, QPainter, QColor
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QVBoxLayout
from PyQt5.QtCore import Qt

# 用于创建一个绘图类
import matplotlib
matplotlib.use("Qt5Agg")  # 声明使用QT5
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import *

import os
import math 

import netCDF4 as nc
import numpy as np

myWin = None

ori_data_matrix = None
fft_data_matrix = None
x = None
y = None
z = None
zero_padded_data_fft_matrix = None
freqs = None
freqs_corr_distance = None

ori_data_width = None 
ori_data_height = None
img_scale = None

current_mouse_x = 0 # 鼠标指向的图片中的x坐标与y坐标
current_mouse_y = 0
zoom_start_distance = 0
zoom_stop_distance = 0

ori_N = 1024 # 原始数据长度
sampling_freq = 1e6 # 采样频率
zero_padded_len = 16384 # 补零后的长度

lut = []

current_file_max_value = 0 # 当前加载的文件中的最大值，用于设置频谱y轴范围 同时用于成像设置
current_file_min_value = 0 # 当前加载的文件中的最小值，用于成像设置
current_lay_freq_index = 0 # 当前显示的图像的层数 Z轴的index

stage_freq_posi_arr = [] # 存储需要叠加显示的数据的坐标
stage_freq_label_arr = []# 存储需要叠加显示的数据的图例

def init_lut():
    '''
    初始化颜色对照表
    此lut可将灰度图映射为伪彩图
    jet
    '''
    global lut
    for i in range(256):
        if(i >= 0 and i <= 31): # 0~31
            r = 0
            g = 0
            b = 128 + 4 * (i - 0)
            tmp_color = QColor(r,g,b)
        elif(i == 32):
            r = 0
            g = 0
            b = 255
            tmp_color = QColor(r,g,b)
        elif(i >= 33 and i <= 95):
            r = 0
            g = 4 + 4 * (i - 33)
            b = 255
            tmp_color = QColor(r,g,b)
        elif(i == 96):
            r = 2
            g = 255
            b = 254
            tmp_color = QColor(r,g,b)
        elif(i >= 97 and i <= 158):
            r = 6 + 4 * (i - 97)
            g = 255
            b = 250 - 4 * (i - 97)
            tmp_color = QColor(r,g,b)
        elif(i == 159):
            r = 254
            g = 255
            b = 1
            tmp_color = QColor(r,g,b)
        elif(i >= 160 and i <= 223):
            r = 255
            g = 252 - 4 * (i - 160) 
            b = 0
            tmp_color = QColor(r,g,b)
        elif(i >= 224 and i <= 255):
            r = 252 - 4 * (i - 224)
            g = 0
            b = 0
            tmp_color = QColor(r,g,b)
        lut.append(tmp_color)
            
def read_data(path):
    dataset = nc.Dataset(path)

    all_vars = dataset.variables.items()

    print("==============================================")
    print("variabbles 数量", len(all_vars))
    print("==============================================")

    all_vars_info = dataset.variables.items()
    print("==============================================")
    print("variables 类型", type(all_vars_info))
    print("==============================================")
    
    all_vars_info = list(all_vars_info)
    print("==============================================")
    for info_ele in all_vars_info:
        print(info_ele)
        print("")
    print("==============================================")


    xx = dataset.variables["xx"]
    yy = dataset.variables["yy"]
    datas = dataset.variables["data"]

    ori_data_array = np.array(datas[:])

    print("len of xx: ", len(xx))
    print("len of yy: ", len(yy))
    print("len of datas: ",len(datas) )

    print(ori_data_array.shape)
    global ori_data_matrix 
    global x, y, z
    x = len(xx)
    y = len(yy)
    z = len(datas)
    print(x)
    print(y)
    print(z)
    ori_data_matrix = np.zeros((x, y, z))

    for tmp_x in range(x):
        for tmp_y in range(y):
            ori_data_matrix[tmp_x, tmp_y, :] = ori_data_array[:,tmp_x , tmp_y]

    print(ori_data_array.shape)
 
def ori_datas_fft():
    # 对信号进行傅里叶变换 进行补零
    global zero_padded_data_fft_matrix
    global freqs
    
    zero_padded_data_matrix = np.zeros((x, y, zero_padded_len))
    zero_padded_data_fft_matrix = np.zeros((x, y, zero_padded_len//2))

    for tmp_x in range(x):
        for tmp_y in range(y):
            zero_padded_data_matrix[tmp_x, tmp_y, :] = np.pad(ori_data_matrix[tmp_x, tmp_y, :], (0, zero_padded_len - ori_N), 'constant')
            tmp_fft_sig = np.fft.fft(zero_padded_data_matrix[tmp_x, tmp_y, :])
            zero_padded_data_fft_matrix[tmp_x, tmp_y, :] = np.abs(tmp_fft_sig)[:zero_padded_len//2]
            
    freqs = np.fft.fftfreq(zero_padded_len, d=1/sampling_freq)
    freqs = freqs[:zero_padded_len//2] # 只取一半
    print('傅里叶变换成果')
    
    # !!!!!!!!!debug
    # 确定方向
    # zero_padded_data_fft_matrix[10:21, 30:41, :] = 0
    # zero_padded_data_fft_matrix[40:51, 21:31, :] = 0
    
def calculate_distance_corr_freqs_index():
    # 计算频率对应的距离关系
    B = 56*10**9      #系统带宽
    T = 1024*10**(-6) #扫频时间
    C = 3*10**11         # 光速
    K = B / T         #扫频斜率
    ref_index = 2.16    # 折射率,根据工件材质不同更改

    global freqs_corr_distance

    freqs_corr_distance = freqs * C / 2 / K / ref_index # 频率对应的距离 #折射率1.2 根据研究对象改变
    oneMM_corr_fre = 2 * K * ref_index / C
    oneMM_corr_fre_index = round(oneMM_corr_fre /(freqs[1] - freqs[0])) 
    print("1mm对应的频点数量 " + str(oneMM_corr_fre_index))
 

 
class MyMainWindow(QMainWindow,Ui_MainWindow): #这里也要记得改
    def __init__(self,parent =None):
        super(MyMainWindow,self).__init__(parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
 
        self.ui.pushButton_select_file.clicked.connect(self.select_file)
        self.ui.pushButton_select_lay.clicked.connect(self.select_lay)
        self.ui.pushButton_move_up.clicked.connect(self.img_up_move)
        self.ui.pushButton_move_down.clicked.connect(self.img_down_move)
        self.ui.pushButton_freq_zoom.clicked.connect(self.freq_zoom)
        
        
        self.ui.label_display_img.setMouseTracking(True)
        self.ui.label_display_img.mouseMoveEvent = self.mouse_move_event
        self.ui.label_display_img.mousePressEvent = self.mousePressEvent
        self.init_freqs_plot_canvas()   
        
    def select_file(self):
        # 初始化一些参数
        global ori_data_width, ori_data_height, img_scale
        ori_data_width = None
        ori_data_height = None
        img_scale = None
        global current_mouse_x, current_mouse_y, zoom_start_distance, zoom_stop_distance
        current_mouse_x = 0
        current_mouse_y = 0
        zoom_start_distance = 0
        zoom_stop_distance = 0
        global current_file_max_value, current_file_min_value, current_lay_freq_index
        current_file_max_value = 0
        current_file_min_value = 0
        current_lay_freq_index = 0
        
        
        # 读取与处理文件
        options = QFileDialog.Options()
        file, _ = QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*);;Text Files (*.txt)", options=options)
        print("select file : " + file)
        read_data(file)
        ori_datas_fft()
        calculate_distance_corr_freqs_index()
        print(zero_padded_data_fft_matrix.shape)
        height, width, _ = zero_padded_data_fft_matrix.shape
        print(height)
        print(width)

        current_file_max_value = np.max(zero_padded_data_fft_matrix)
        
        ori_data_width = width
        ori_data_height = height
        
        
        # 获取这个文件中的最大值与最小值
        current_file_max_value = np.max(zero_padded_data_fft_matrix)
        current_file_min_value = np.min(zero_padded_data_fft_matrix)
        
        # 将本层的数据归一化，按照本层的最大值与最小值来
        tmp_max = np.max(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
        tmp_min = np.min(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
        normalized_data = ((zero_padded_data_fft_matrix[:,:,current_lay_freq_index] - tmp_min) / (tmp_max - tmp_min)) * 255
        normalized_data = normalized_data.astype(np.uint8)
        
        image = QImage(ori_data_width, ori_data_height, QImage.Format_RGB32)
        for w in range(image.width()):
            for h in range(image.height()):
                image.setPixel(w, h, lut[normalized_data[h,w]].rgb())
        pixmap = QPixmap.fromImage(image)
        
        # 缩放 将图片缩放到合适的大小
        tmp_label_width = self.ui.label_display_img.width()
        tmp_label_heigh = self.ui.label_display_img.height()
        
        # 定义缩放倍率 = 缩放后 / 原始图片
        # 缩放后 = 原始图片 * 缩放倍率
        # 原始图 = 缩放后 / 缩放倍率
        tmp_x_sacle = tmp_label_width / width
        tmp_y_scale = tmp_label_heigh / height
        
        
        global img_scale_width, img_scale_height
        img_scale = min(tmp_x_sacle, tmp_y_scale)
        img_scale = math.floor(img_scale)
        img_scale_width = width * img_scale
        img_scale_height = height * img_scale
        
        scaled_pixmap = pixmap.scaled(img_scale_width, img_scale_height, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    
        self.ui.label_display_img.setPixmap(scaled_pixmap)
        self.ui.label_display_img.setAlignment(Qt.AlignLeft | Qt.AlignTop) # 使图片显示于左上角
        
    def select_lay(self):
        global current_lay_freq_index
        
        self.ui.label_display_img.clear()
        
        current_distance =  self.ui.spinBox_select_distance.value()
        # 计算距离对应的频率
        differences = np.abs(freqs_corr_distance - current_distance)
        current_lay_freq_index = np.argmin(differences)
        
        # 如果设置了频谱的显示范围，就按照频谱的显示范围确定灰度的最大值与最小值
        # 如果没有设置频谱显示范围，就按照本层数据决定灰度的最大值与最小值
        tmp_max = 0
        tmp_min = 0
        if zoom_start_distance != 0 or zoom_stop_distance != 0:
            differences = np.abs(freqs_corr_distance - zoom_start_distance)
            freq_start_index = np.argmin(differences)
            differences = np.abs(freqs_corr_distance - zoom_stop_distance)
            freq_stop_index = np.argmin(differences)
            tmp_max = np.max(zero_padded_data_fft_matrix[:, :, freq_start_index:freq_stop_index])
            tmp_min = np.min(zero_padded_data_fft_matrix[:, :, freq_start_index:freq_stop_index])
        else:
            tmp_max = np.max(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
            tmp_min = np.min(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
            
        normalized_data = ((zero_padded_data_fft_matrix[:,:,current_lay_freq_index] - tmp_min) / (tmp_max - tmp_min)) * 255
        normalized_data = normalized_data.astype(np.uint8)
        # image = QImage(normalized_data.data, ori_data_width, ori_data_height, ori_data_width, QImage.Format_RGB32)
        image = QImage(ori_data_width, ori_data_height, QImage.Format_RGB32)
        for w in range(image.width()):
            for h in range(image.height()):
                image.setPixel(w, h, lut[normalized_data[h,w]].rgb())
        
        pixmap = QPixmap.fromImage(image)
        scaled_pixmap = pixmap.scaled(img_scale_width, img_scale_height, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    
        self.ui.label_display_img.setPixmap(scaled_pixmap)
        self.ui.label_display_img.setAlignment(Qt.AlignLeft | Qt.AlignTop) # 使图片显示于左上角
        
    def mouse_move_event(self, event):
        if ori_data_matrix is None:
            return
        """捕获鼠标移动事件并显示坐标"""
        # 获取鼠标相对于 QLabel 的位置
        position = event.pos()

        # 获取鼠标在图片中的坐标
        x = position.x()
        y = position.y()

        tmp_x = math.floor(x / img_scale) 
        tmp_y = math.floor(y / img_scale) 
        # 检查鼠标是否在图片范围内
        if 0 <= tmp_x < ori_data_width and 0 <= tmp_y < ori_data_height:
            # print(f"鼠标指向的图片坐标: ({x}, {y})")
            self.ui.label_pointer_x.setText(str(tmp_x))
            self.ui.label_pointer_y.setText(str(tmp_y))
            self.ui.label_pointer_current_value.setText(str(zero_padded_data_fft_matrix[tmp_y, tmp_x, 0]))
            
            global current_mouse_x, current_mouse_y
            current_mouse_x = tmp_x
            current_mouse_y = tmp_y
            
            self.plot_freq()
        else:
            self.ui.label_pointer_x.setText("out")
            self.ui.label_pointer_y.setText("out")
            self.ui.label_pointer_current_value.setText("out")
    
    def mousePressEvent(self, event):
        global stage_freq_label_arr, stage_freq_posi_arr
        mouse_position = event.pos()
        tmp_x = mouse_position.x()
        tmp_y = mouse_position.y()
        tmp_x = math.floor(tmp_x / img_scale) 
        tmp_y = math.floor(tmp_y / img_scale) 
        if 0 <= tmp_x < ori_data_width and 0 <= tmp_y < ori_data_height:
            text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter label info:')
            if ok and text:
                tmp_posi = [tmp_x, tmp_y]
                stage_freq_posi_arr.append(tmp_posi)
                stage_freq_label_arr.append(text)
            elif not text:
                return
            
            
    
    def init_freqs_plot_canvas(self):
        tmp_w = self.ui.widget_display_freq.width()
        tmp_h = self.ui.widget_display_freq.height()
        self.F = MyFigure(width = tmp_w, height = tmp_h, dpi=100)
        self.F.fig.suptitle("Amplitude-Distance")
        self.F.axes.grid(True)
        
        self.gridlayout = QGridLayout(self.ui.widget_display_freq)  # 继承容器groupBox
        self.gridlayout.addWidget(self.F,0,1)
    
    def plot_freq(self):
        # print(str(tmp_x) + " " + str(tmp_y))
        self.F.axes.clear()
        
        # 优先显示需要显示的位置的数据
        if len(stage_freq_label_arr) != 0 :
            if zoom_start_distance !=0 or zoom_stop_distance != 0:   
                differences = np.abs(freqs_corr_distance - zoom_start_distance)
                start_index = np.argmin(differences)
                
                differences = np.abs(freqs_corr_distance - zoom_stop_distance)
                stop_index = np.argmin(differences)
                for i in range(len(stage_freq_label_arr)):
                    tmp_x = stage_freq_posi_arr[i][0]
                    tmp_y = stage_freq_posi_arr[i][1]
                    
                    self.F.axes.plot(freqs_corr_distance[start_index:stop_index], zero_padded_data_fft_matrix[tmp_y, tmp_x, start_index:stop_index], label=stage_freq_label_arr[i])
                    self.F.axes.set_ylim(0, current_file_max_value)
                self.F.axes.legend()
            else:
                for i in range(len(stage_freq_label_arr)):
                    tmp_x = stage_freq_posi_arr[i][0]
                    tmp_y = stage_freq_posi_arr[i][1]
                    
                    self.F.axes.plot(freqs_corr_distance, zero_padded_data_fft_matrix[tmp_y, tmp_x, :], label=stage_freq_label_arr[i])
                    self.F.axes.set_ylim(0, current_file_max_value)
                self.F.axes.legend()
                
        
        if zoom_start_distance !=0 or zoom_stop_distance != 0:            
            differences = np.abs(freqs_corr_distance - zoom_start_distance)
            start_index = np.argmin(differences)
            
            differences = np.abs(freqs_corr_distance - zoom_stop_distance)
            stop_index = np.argmin(differences)
            self.F.axes.plot(freqs_corr_distance[start_index:stop_index], zero_padded_data_fft_matrix[current_mouse_y, current_mouse_x, start_index:stop_index], 'r-')
            self.F.axes.set_ylim(0, current_file_max_value)
        
        else:
            self.F.axes.plot(freqs_corr_distance, zero_padded_data_fft_matrix[current_mouse_y, current_mouse_x, :], 'r-')
            self.F.axes.set_ylim(0, current_file_max_value)

        self.F.draw()

    def img_up_move(self):
        global current_lay_freq_index
        if current_lay_freq_index > 0:
            current_lay_freq_index = current_lay_freq_index - 1
        
        tmp_max = 0
        tmp_min = 0
        if zoom_start_distance != 0 or zoom_stop_distance != 0:
            differences = np.abs(freqs_corr_distance - zoom_start_distance)
            freq_start_index = np.argmin(differences)
            differences = np.abs(freqs_corr_distance - zoom_stop_distance)
            freq_stop_index = np.argmin(differences)
            tmp_max = np.max(zero_padded_data_fft_matrix[:, :, freq_start_index:freq_stop_index])
            tmp_min = np.min(zero_padded_data_fft_matrix[:, :, freq_start_index:freq_stop_index])
        else:
            tmp_max = np.max(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
            tmp_min = np.min(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
            
        normalized_data = ((zero_padded_data_fft_matrix[:,:,current_lay_freq_index] - tmp_min) / (tmp_max - tmp_min)) * 255
        normalized_data = normalized_data.astype(np.uint8)
        # image = QImage(normalized_data.data, ori_data_width, ori_data_height, ori_data_width, QImage.Format_RGB32)
        image = QImage(ori_data_width, ori_data_height, QImage.Format_RGB32)
        for w in range(image.width()):
            for h in range(image.height()):
                image.setPixel(w, h, lut[normalized_data[h,w]].rgb())
        
        pixmap = QPixmap.fromImage(image)
        scaled_pixmap = pixmap.scaled(img_scale_width, img_scale_height, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    
        self.ui.label_display_img.setPixmap(scaled_pixmap)
        self.ui.label_display_img.setAlignment(Qt.AlignLeft | Qt.AlignTop) # 使图片显示于左上角
        
    
    def img_down_move(self):
        global current_lay_freq_index
        if current_lay_freq_index < len(freqs_corr_distance) - 1:
            current_lay_freq_index = current_lay_freq_index + 1
        
        tmp_max = 0
        tmp_min = 0
        if zoom_start_distance != 0 or zoom_stop_distance != 0:
            differences = np.abs(freqs_corr_distance - zoom_start_distance)
            freq_start_index = np.argmin(differences)
            differences = np.abs(freqs_corr_distance - zoom_stop_distance)
            freq_stop_index = np.argmin(differences)
            
            print(zoom_start_distance)
            print(zoom_stop_distance)
            print(freq_start_index)
            print(freq_stop_index)
            
            tmp_max = np.max(zero_padded_data_fft_matrix[:, :, freq_start_index:freq_stop_index])
            tmp_min = np.min(zero_padded_data_fft_matrix[:, :, freq_start_index:freq_stop_index])
        else:
            tmp_max = np.max(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
            tmp_min = np.min(zero_padded_data_fft_matrix[:, :, current_lay_freq_index])
            
        normalized_data = ((zero_padded_data_fft_matrix[:,:,current_lay_freq_index] - tmp_min) / (tmp_max - tmp_min)) * 255
        normalized_data = normalized_data.astype(np.uint8)
        # image = QImage(normalized_data.data, ori_data_width, ori_data_height, ori_data_width, QImage.Format_RGB32)
        image = QImage(ori_data_width, ori_data_height, QImage.Format_RGB32)
        for w in range(image.width()):
            for h in range(image.height()):
                image.setPixel(w, h, lut[normalized_data[h,w]].rgb())
        
        pixmap = QPixmap.fromImage(image)
        scaled_pixmap = pixmap.scaled(img_scale_width, img_scale_height, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    
        self.ui.label_display_img.setPixmap(scaled_pixmap)
        self.ui.label_display_img.setAlignment(Qt.AlignLeft | Qt.AlignTop) # 使图片显示于左上角
    
    def freq_zoom(self):
        global zoom_start_distance, zoom_stop_distance
        zoom_start_distance = self.ui.spinBox_zoom_start_distance.value()
        zoom_stop_distance = self.ui.spinBox_zoom_stop_distance.value()
        
        differences = np.abs(freqs_corr_distance - zoom_start_distance)
        start_index = np.argmin(differences)
        
        differences = np.abs(freqs_corr_distance - zoom_stop_distance)
        stop_index = np.argmin(differences)
        
        self.F.axes.clear()
        self.F.axes.plot(freqs_corr_distance[start_index : stop_index], \
            zero_padded_data_fft_matrix[current_mouse_y, current_mouse_x, start_index:stop_index], 'r-')
        self.F.draw()
        
class MyFigure(FigureCanvas):
    def __init__(self,width=5, height=4, dpi=100):
        #第一步：创建一个创建Figure
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        #第二步：在父类中激活Figure窗口
        super(MyFigure,self).__init__(self.fig) #此句必不可少，否则不能显示图形
        #第三步：创建一个子图，用于绘制图形用，111表示子图编号，如matlab的subplot(1,1,1)
        self.axes = self.fig.add_subplot(111)
        
 
if __name__ == "__main__":
    app = QApplication(sys.argv)
    init_lut() # 初始化伪彩图颜色对照表
    myWin = MyMainWindow()

    
    
    myWin.show()
    sys.exit(app.exec_())    