__author__ = 'sean'
import Tkinter as tk
import multiprocessing as mp
from Tkinter import Tk, Button, Text, Checkbutton, IntVar, N, S, E, W
from RecreateData import *
from GUIMisc import *
from RecreateDataGUISettings import *
import sys, os, webbrowser, re, cPickle, math, time, __future__

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def wrapped_partial(func, *args, **kwargs):
    partial_func = functools.partial(func, *args, **kwargs)
    functools.update_wrapper(partial_func, func)
    return partial_func

def doc_call_wrapper(name, doc_dir):
    if webbrowser._tryorder[0] == 'MacOSX' and "chrome" not in webbrowser._tryorder and (webbrowser._tryorder)>1: # This is a chrome bug work-around to make sure we aren't using Chrome on Os High Serria, delete the if statement at a later date.
        return lambda a=0: webbrowser.get(webbrowser._tryorder[1]).open("file://" + doc_dir + "/corvids_readme.html#" + name, new=2)
    return lambda a=0: webbrowser.get().open("file://" + doc_dir + "/corvids_readme.html#" + name, new=2)


class PossValsText(tk.Text):

    '''A text widget that accepts a 'min_value' and 'max_value' option'''

    def __init__(self, parent, *args, **kwargs):
        try:
            self._min_val = kwargs.pop("min_value")
            self._max_val = kwargs.pop("max_value")
        except KeyError:
            raise SyntaxError

        tk.Text.__init__(self, parent, *args, **kwargs)

        # if the variable has data in it, use it to initialize
        # the widget
        if self._min_val is not None and self._max_val is not None:
            var_current = [str(n) for n in xrange(self._min_val.get(), self._max_val.get()+1)]
            self.insert("1.0", ', '.join(var_current))

        # this defines an internal proxy which generates a
        # virtual event whenever text is inserted or deleted
        self.tk.eval('''
            proc widget_proxy {widget widget_command args} {

                # call the real tk widget command with the real args
                set result [uplevel [linsert $args 0 $widget_command]]

                # if the contents changed, generate an event we can bind to
                if {([lindex $args 0] in {insert replace delete})} {
                    event generate $widget <<Change>> -when tail
                }
                # return the result from the real widget command
                return $result
            }
            ''')

        # this replaces the underlying widget with the proxy
        self.tk.eval('''
            rename {widget} _{widget}
            interp alias {{}} ::{widget} {{}} widget_proxy {widget} _{widget}
        '''.format(widget=str(self)))

        # # set up a binding to update the variable whenever
        # # the widget changes
        # self.bind("<<Change>>", self._on_widget_change)

        # set up a trace to update the text widget when the
        # variable changes
        if self._min_val is not None:
            self._min_val.trace("wu", self._on_var_change)
        if self._max_val is not None:
            self._max_val.trace("wu", self._on_var_change)

    def _on_var_change(self, *args):
        '''Change the text widget when the associated textvariable changes'''

        # only change the widget if something actually
        # changed, otherwise we'll get into an endless
        # loop
        try:
            self._min_val.get()
            self._max_val.get()
        except ValueError:
            return
        current_poss = re.split('[,\s]+', self.get("1.0", "end-1c"))
        var_current = [str(n) for n in xrange(self._min_val.get(), self._max_val.get()+1)]
        if current_poss != var_current:
            self.delete("1.0", "end")
            self.insert("1.0", ', '.join(var_current))

    def getPossVals(self):
        current_poss = re.split('[,\s]+', self.get("1.0", "end-1c"))
        poss_vals = [int(math.floor(eval(poss))) for poss in current_poss if
                len(poss)
                and (poss=='0' or (poss if poss.find('..') > -1
                                     else poss.lstrip('-+').rstrip('0').rstrip('.')).isdigit())
                and math.floor(eval(poss))==eval(poss)]
        # poss_vals = [int(poss) for poss in current_poss]
        poss_vals.sort()
        return poss_vals

    # def _on_widget_change(self, event=None):
    #     '''Change the variable when the widget changes'''
    #     if self._textvariable is not None:
    #         self._textvariable.set(self.get("1.0", "end-1c"))


class StdoutRedirector(object):

    def __init__(self, text_widget):
        self.set = False
        self.text_space = text_widget
        self.text = ""

    def write(self,string):
        self.text += string
        self.text_space.insert('end', string)
        self.text_space.see('end')

    def flush(self):
        pass

class CoreGUI(object):

    def __init__(self,parent, doc_dir):
        self.doc_dir = doc_dir
        self.parent = parent
        self.InitUI()
        self.show_poss_vals = False
        self.show_forced_vals = False
        self.use_SD = False
        self.out = StdoutRedirector(self.text_box)
        sys.stdout = self.out

    def main(self):
        '''
        The body of the execution.  It is split up so that the GUI freezes as little as possible, but it will still
        freeze during several portions (some potentially quite long)
        :return:
        '''
        self.start_button.config(text='Running...')
        self.start_button.config(command=None)
        self.graph_button.grid_remove()
        self.start_button.grid_remove()
        self.start_button.grid(row=15, columnspan=3, column=0)
        debug = bool(self.debug.get())
        if debug:
            print "Starting Data Recreation..."
        var = eval(compile(self.variance.get(), '<string>', 'eval', __future__.division.compiler_flag))
        var_precision = abs(eval(compile(self.variance_precision.get(), '<string>', 'eval', __future__.division.compiler_flag)))
        if self.use_SD:
            min_var = max(0, (var - var_precision))**2
            max_var = (var + var_precision)**2
            var = (min_var + max_var)/2
            var_precision = (max_var - min_var)/2


        self.rd = RecreateData(self.min_score.get(), self.max_score.get(), self.num_subjects.get(),
                               eval(compile(self.mean.get(), '<string>', 'eval', __future__.division.compiler_flag)),
                               var,
                               debug=debug,
                               mean_precision=eval(compile(self.mean_precision.get(), '<string>', 'eval', __future__.division.compiler_flag)),
                               variance_precision=var_precision
                               )
        # self.rd.recreateData(check_val=self.getForcedVals(), poss_vals=self.poss_vals.getPossVals(),
        #                             multiprocess=True, find_first=(not bool(self.find_all.get())))
        self.parent.after(5, self.main_2)

    def main_2(self):
        debug = bool(self.debug.get())
        mean_var_pairs = self.rd._recreateData_piece_1(check_val=self.getForcedVals(), poss_vals=self.poss_vals.getPossVals(),
                                    multiprocess=True, find_first=(bool(self.find_first.get())))
        partial_3 = wrapped_partial(self.main_3, mean_var_pairs)
        self.parent.after(5, partial_3)

    def main_3(self, mean_var_pairs):
        debug = bool(self.debug.get())
        solution_spaces = self.rd._recreateData_piece_2(mean_var_pairs, check_val=self.getForcedVals(), poss_vals=self.poss_vals.getPossVals(),
                                    multiprocess=True, find_first=(bool(self.find_first.get())))
        partial_4 = wrapped_partial(self.main_4, solution_spaces)
        self.parent.after(5, partial_4)

    def main_4(self, solution_spaces):
        debug = bool(self.debug.get())
        if (bool(self.find_first.get())):
            self.rd._findFirst_piece_1(solution_spaces, check_val=self.getForcedVals(), poss_vals=self.poss_vals.getPossVals(),
                                    multiprocess=True, find_first=(bool(self.find_first.get())))
            self.parent.after(5, self.main_6)
        else:
            init_bases, init_base_vecs, param_tuples = self.rd._findAll_piece_1_multi_proc(solution_spaces, check_val=self.getForcedVals(), poss_vals=self.poss_vals.getPossVals(),
                                    multiprocess=True, find_first=(bool(self.find_first.get())))
            partial_5 = wrapped_partial(self.main_5, init_bases, init_base_vecs, param_tuples)
            self.parent.after(5, partial_5)

    def main_5(self, init_bases, init_base_vecs, param_tuples):
        debug = bool(self.debug.get())
        self.rd._findAll_piece_2_multi_proc(init_bases, init_base_vecs, param_tuples, check_val=self.getForcedVals(), poss_vals=self.poss_vals.getPossVals(),
                                    multiprocess=True, find_first=(bool(self.find_first.get())))
        self.parent.after(5, self.main_6)

    def main_6(self):
        debug = bool(self.debug.get())
        self.start_button.config(text='Start')
        self.start_button.config(command=self.main)
        self.graph_button.grid()
        self.start_button.grid_remove()
        self.start_button.grid(row=15, columnspan=1, column=0)
        if len(self.rd.sols)>0:
            self.rd.getDataSimple()
            if debug:
                print str(sum([len(x) for x in self.rd.simpleData.itervalues()])) + " unique solutions found."
                index = 0
                for params in self.rd.simpleData:
                    if index > 100:
                        break
                    print "At mean, variance", params, ":"
                    for simpleSol in self.rd.simpleData[params]:
                        if index > 100:
                            break
                        index += 1
                        print simpleSol
        else:
            if debug:
                print str(len(self.rd.getDataSimple())) + " unique solutions found."

    def graph(self):
        self.rd.graphData()

    def _showOrHidePossVals(self):
        self.show_poss_vals = not self.show_poss_vals
        if self.show_poss_vals:
            self.poss_vals.grid()
        else:
            self.poss_vals.grid_remove()

    def _showOrHideForcedVals(self):
        self.show_forced_vals = not self.show_forced_vals
        if self.show_forced_vals:
            self.forced_vals.grid()
        else:
            self.forced_vals.grid_remove()

    def _switch_SD_Var(self):
        self.use_SD = not self.use_SD
        if self.use_SD:
            self.SD_Button.config(text='Switch to Variance Instead of SD')
            self.variance_label.config(text='Standard Deviation')
            self.variance.set(0.0)
            self.variance_precision_label.config(text='SD Precision')
            self.variance_precision.set(0.0)
        else:
            self.SD_Button.config(text='Switch to SD Instead of Variance')
            self.variance_label.config(text='Variance')
            self.variance.set(0.0)
            self.variance_precision_label.config(text='Variance Precision')
            self.variance_precision.set(0.0)

    def getForcedVals(self):
        current_forced = re.split('[,\s]+', self.forced_vals.get("1.0", "end-1c"))
        forced_vals = [int(math.floor(eval(forced))) for forced in current_forced if
                len(forced)
                and (forced=='0' or (forced if forced.find('..') > -1
                                     else forced.lstrip('-+').rstrip('0').rstrip('.')).isdigit())
                and math.floor(eval(forced))==eval(forced)]
        forced_vals.sort()
        return forced_vals

    def clearInputFields(self):
        self.variance.set(0.0)
        self.mean.set(0.0)
        self.min_score.set(0)
        self.max_score.set(0)
        self.variance_precision.set(0.0)
        self.mean_precision.set(0.0)
        self.num_subjects.set(0)

    def clearOutputField(self):
        self.text_box.delete('1.0', END)

    def clearAll(self):
        self.clearInputFields()
        self.clearOutputField()

    def InitUI(self):

        from PIL import Image
        from PIL import ImageTk

        self.parent.lift()
        # image = Image.open("Corvid.png")
        # help_image = Image.open("Question.png")
        # photo = ImageTk.PhotoImage(image)
        # help_photo = ImageTk.PhotoImage(help_image)
        # img = PhotoImage(file="Corvid.png")
        # self.parent.tk.call('wm', 'iconphoto', self.parent._w, photo)
        self.parent.wm_title(string="CORVIDS")

        self.window = MenuItems(self.parent)
        self.window.addMenu("clearMenu", "Clear")
        self.window.addTopLevelMenuItem(self.window.fileMenu, "Load Settings", self.myLoadFunc(self.window), hotKey="<Control-o>")
        self.window.addTopLevelMenuItem(self.window.fileMenu, "Save Settings", self.mySaveFunc(self.window), hotKey="<Control-s>")
        self.window.addTopLevelMenuItem(self.window.fileMenu, "Load Model", self.myLoadModelFunc(self.window), hotKey="<Control-m>")
        self.window.addTopLevelMenuItem(self.window.fileMenu, "Save Model", self.mySaveModelFunc(self.window), hotKey="<Control-M>")
        self.window.addTopLevelMenuItem(self.window.fileMenu, "Save Data", self.mySaveDataFunc(self.window), hotKey="<Control-d>")
        self.window.addTopLevelMenuItem(self.window.clearMenu, "Clear Input Fields", self.clearInputFields)
        self.window.addTopLevelMenuItem(self.window.clearMenu, "Clear Output Field", self.clearOutputField)
        self.window.addTopLevelMenuItem(self.window.clearMenu, "Clear All", self.clearAll)

        self.window.addQuit()

        description_col = 0
        field_col = 2
        help_col = 1

        self.debug = IntVar()
        self.debug.set(1)
        Checkbutton(self.parent, text="Show Progress", variable=self.debug).grid(column=description_col, row=0, sticky=W)
        self.find_first = IntVar()
        self.find_first.set(0)
        Checkbutton(self.parent, text="Stop After First Solution", variable=self.find_first).grid(column=field_col, row=0, sticky=E)

        Label(self.parent, text="Minimum Value", anchor='e').grid(row=1, column=description_col, sticky=N+S+E+W)
        self.min_score = IntVar()
        self.min_entry = Entry(self.parent, textvariable=self.min_score)
        self.min_entry.grid(row=1, column=field_col, sticky=E+W, padx=15)
        min_value_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("minimum-value", self.doc_dir)) #width=14, height=14)#, text="?")
        # min_value_help.config(image=help_photo)
        # min_value_help.image = help_photo
        min_value_help.grid(row=1, column=help_col, sticky=E+W)


        Label(self.parent, text="Maximum Value", anchor='e').grid(row=2, column=description_col, sticky=N+S+E+W)
        self.max_score = IntVar()
        self.max_entry = Entry(self.parent, textvariable=self.max_score)
        self.max_entry.grid(row=2, column=field_col, sticky=E+W, padx=15)
        max_value_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("maximum-value", self.doc_dir)) #width=14, height=14)#, text="?")
        # max_value_help.config(image=help_photo)
        # max_value_help.image = help_photo
        max_value_help.grid(row=2, column=help_col, sticky=E+W)

        Label(self.parent, text="Mean", anchor='e').grid(row=3, column=description_col, sticky=N+S+E+W)
        self.mean = StringVar()
        self.mean.set(0.0)
        self.mean_entry = Entry(self.parent, textvariable=self.mean)
        self.mean_entry.grid(row=3, column=field_col, sticky=E+W, padx=15)
        mean_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("mean", self.doc_dir)) #width=14, height=14)#, text="?")
        # mean_help.config(image=help_photo)
        # mean_help.image = help_photo
        mean_help.grid(row=3, column=help_col, sticky=E+W)

        Label(self.parent, text="Mean Precision", anchor='e').grid(row=4, column=description_col, sticky=N+S+E+W)
        self.mean_precision = StringVar()
        self.mean_precision.set(0.0)
        self.mean_precision_entry = Entry(self.parent, textvariable=self.mean_precision)
        self.mean_precision_entry.grid(row=4, column=field_col, sticky=E+W, padx=15)
        mean_precision_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("mean-precision", self.doc_dir)) #width=14, height=14)#, text="?")
        # mean_precision_help.config(image=help_photo)
        # mean_precision_help.image = help_photo
        mean_precision_help.grid(row=4, column=help_col, sticky=E+W)

        self.SD_Button = Button(self.parent, text="Switch to SD Instead of Variance", command=self._switch_SD_Var)
        self.SD_Button.grid(row=5, columnspan=3, column=description_col, sticky=E+W)

        self.variance_label = Label(self.parent, text="Variance", anchor='e')
        self.variance_label.grid(row=6, column=description_col, sticky=N+S+E+W)
        self.variance = StringVar()
        self.variance.set(0.0)
        self.variance_entry = Entry(self.parent, textvariable=self.variance)
        self.variance_entry.grid(row=6, column=field_col, sticky=E+W, padx=15)
        variance_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("variance", self.doc_dir)) #width=14, height=14)#, text="?")
        # variance_help.config(image=help_photo)
        # variance_help.image = help_photo
        variance_help.grid(row=6, column=help_col, sticky=E+W)

        self.variance_precision_label = Label(self.parent, text="Variance Precision", anchor='e')
        self.variance_precision_label.grid(row=7, column=description_col, sticky=N+S+E+W)
        self.variance_precision = StringVar()
        self.variance_precision.set(0.0)
        self.variance_precision_entry = Entry(self.parent, textvariable=self.variance_precision)
        self.variance_precision_entry.grid(row=7, column=field_col, sticky=E+W, padx=15)
        variance_precision_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("variance-precision", self.doc_dir)) #width=14, height=14)#, text="?")
        # variance_precision_help.config(image=help_photo)
        # variance_precision_help.image = help_photo
        variance_precision_help.grid(row=7, column=help_col, sticky=E+W)

        Label(self.parent, text="Number of Subjects", anchor='e').grid(row=8, column=description_col, sticky=N+S+E+W)
        self.num_subjects = IntVar()
        self.num_subjects_entry = Entry(self.parent, textvariable=self.num_subjects)
        self.num_subjects_entry.grid(row=8, column=field_col, sticky=E+W, padx=15)
        num_subjects_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("number-of-subjects", self.doc_dir)) #width=14, height=14)#, text="?")
        # num_subjects_help.config(image=help_photo)
        # num_subjects_help.image = help_photo
        num_subjects_help.grid(row=8, column=help_col, sticky=E+W)



        Label(self.parent, text="Possible Values", anchor='e').grid(row=9, column=description_col, sticky=N+S+E+W)
        # Button(self.parent, text="Possible Values", command=self._showOrHidePossVals).grid(row=8, columnspan=2, column=description_col, sticky=E+W)
        Button(self.parent, text="Show/Hide", command=self._showOrHidePossVals).grid(row=9, column=field_col, sticky=E+W, padx=15)
        self.poss_vals = PossValsText(self.parent, min_value=self.min_score, max_value=self.max_score, wrap='word', height = 11, width=50)
        self.poss_vals.grid(row=10, columnspan=3, column=description_col, sticky=N+S+E+W)
        self.poss_vals.grid_remove()
        # poss_vals_frame = Frame(self.parent, height=10, width=10)
        # poss_vals_frame.grid(row=10, column=help_col, sticky=E+W)
        # poss_vals_help = Button(poss_vals_frame, text="?", font=('TkDefaultFont', 8)) #width=14, height=16)#, text="?")
        # poss_vals_help.pack()
        poss_vals_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("possible-values", self.doc_dir)) #width=14, height=16)#, text="?")
        # poss_vals_help.config(image=help_photo)
        # poss_vals_help.image = help_photo
        poss_vals_help.grid(row=9, column=help_col, sticky=E+W)

        Label(self.parent, text="Forced Values", anchor='e').grid(row=11, column=description_col, sticky=N+S+E+W)
        Button(self.parent, text="Show/Hide", command=self._showOrHideForcedVals).grid(row=11, column=field_col, sticky=E+W, padx=15)
        # Button(self.parent, text="Forced Values", command=self._showOrHideForcedVals).grid(row=10, columnspan=2, column=description_col, sticky=E+W)
        self.forced_vals = Text(self.parent, wrap='word', height = 12, width=50)
        self.forced_vals.grid(row=12, columnspan=3, column=description_col, sticky=N+S+E+W)
        self.forced_vals.grid_remove()
        # forced_vals_frame = Frame(self.parent, height=10, width=10)
        # forced_vals_frame.grid(row=10, column=help_col, sticky=E+W)
        forced_vals_help = Button(self.parent, text="?", font=('TkDefaultFont', 8), command=doc_call_wrapper("forced-values", self.doc_dir)) #width=14, height=16)#, text="?")
        # forced_vals_help.config(image=help_photo)
        # forced_vals_help.image = help_photo
        forced_vals_help.grid(row=11, column=help_col, sticky=E+W)

        Label(self.parent, text="Output:").grid(row=13, columnspan=3, column=description_col, sticky=E+W)
        # self.outputFrame = tk.Frame(self.parent)
        self.parent.columnconfigure(0,weight=1)
        self.parent.rowconfigure(14, weight=1)
        self.text_box = Text(self.parent, wrap='word', height = 11, width=50)
        self.text_box.grid(row=14, columnspan=3, column=0, sticky=N+S+E+W)

        self.start_button = Button(self.parent, text="Start", command=self.main)
        self.start_button.grid(row=15, columnspan=3, column=0)

        self.graph_button = Button(self.parent, text="Graph", command=self.graph)
        self.graph_button.grid(row=15, column=field_col)
        self.graph_button.grid_remove()


    def setSettings(self, settings_file_full_path):
        '''
        Takes a path and a filename and sets the simulation to that, informing the user that it has been loaded.
        :param settings_file_full_path: from root directory to filename and extension
        :return:
        '''
        if self.debug.get():
            print "Loading settings from " + settings_file_full_path
        self.settings_file = settings_file_full_path
        input_file = open(self.settings_file, 'rb')
        settings = cPickle.load(input_file)
        assert isinstance(settings, RecreateDataGUISettings)
        input_file.close()
        self.debug.set(settings.debug)
        self.min_score.set(settings.min_score)
        self.max_score.set(settings.max_score)
        self.poss_vals.delete("1.0", "end")
        self.poss_vals.insert("1.0", settings.poss_vals)
        self.mean.set(settings.mean)
        self.mean_precision.set(settings.mean_precision)
        if self.use_SD != settings.use_SD:
            self._switch_SD_Var()
        self.variance.set(settings.variance)
        self.variance_precision.set(settings.variance_precision)
        self.num_subjects.set(settings.num_samples)
        self.forced_vals.delete("1.0", "end")
        self.forced_vals.insert("1.0", settings.check_vals)


    def saveSettings(self, settings_file_full_path):
        if self.debug.get():
            print "Saving settings to " + settings_file_full_path
        self.settings_file = settings_file_full_path
        settings_save_file = open(self.settings_file, 'wb')
        settings = RecreateDataGUISettings(self.debug.get(), self.min_score.get(), self.max_score.get(), self.poss_vals.get("1.0", "end-1c"), self.mean.get(), self.mean_precision.get(), self.variance.get(), self.variance_precision.get(), self.num_subjects.get(), self.forced_vals.get("1.0", "end-1c"), self.use_SD)
        cPickle.dump(settings, settings_save_file)
        settings_save_file.close()

    def saveModel(self, model_file_full_path):
        if self.debug.get():
            print "Saving model to " + model_file_full_path
        self.model_file = model_file_full_path
        model_save_file = open(self.model_file, 'wb')
        cPickle.dump(self.rd, model_save_file)
        model_save_file.close()

    def loadModel(self, model_file_full_path):
        if self.debug.get():
            print "Loading model from " + model_file_full_path

        self.model_file = model_file_full_path
        model_read_file = open(self.model_file, 'rb')
        self.rd = cPickle.load(model_read_file)
        assert isinstance(self.rd, RecreateData)
        model_read_file.close()
        self.graph_button.grid()
        self.start_button.grid_remove()
        self.start_button.grid(row=14, columnspan=1, column=0)

    def saveData(self, data_file_full_path):
        if self.debug.get():
            print "Saving data to " + data_file_full_path
        self.data_file = data_file_full_path
        data_save_file = open(self.data_file, 'wb')
        if len(self.rd.sols) == 0:
            data_save_file.close()
            return
        dataStr = ""
        for param, solutions in self.rd.simpleData.iteritems():
            dataStr += str(param) + ":\n"
            for sol in solutions:
                dataStr += "\t" + str(sol) + "\n"
        data_save_file.write(dataStr.strip())
        data_save_file.close()


    def myLoadFunc_lambda(self,environment):
        lambda_new_window("openerWindow",environment, "Load Settings")
        lambda_make_file_loader(environment,"openerWindow","opener",self.setSettings)

    def myLoadFunc(self,environment):
        return lambda environment=environment: self.myLoadFunc_lambda(environment)

    def mySaveFunc_lambda(self,environment):
        lambda_new_window("openerWindow",environment, "Save Settings")
        lambda_make_file_saver(environment,"openerWindow","opener",self.saveSettings)

    def mySaveFunc(self,environment):
        return lambda environment=environment: self.mySaveFunc_lambda(environment)

    def mySaveDataFunc_lambda(self,environment):
        lambda_new_window("openerWindow",environment, "Save Data")
        lambda_make_file_saver(environment,"openerWindow","opener",self.saveData)

    def mySaveDataFunc(self,environment):
        return lambda environment=environment: self.mySaveDataFunc_lambda(environment)

    def mySaveModelFunc_lambda(self,environment):
        lambda_new_window("openerWindow",environment, "Save Model")
        lambda_make_file_saver(environment,"openerWindow","opener",self.saveModel)

    def mySaveModelFunc(self,environment):
        return lambda environment=environment: self.mySaveModelFunc_lambda(environment)

    def myLoadModelFunc_lambda(self,environment):
        lambda_new_window("openerWindow",environment, "Load Model")
        lambda_make_file_loader(environment,"openerWindow","opener",self.loadModel)

    def myLoadModelFunc(self,environment):
        return lambda environment=environment: self.myLoadModelFunc_lambda(environment)
def root_closer(root, directory_name=None):
    if directory_name:
        shutil.rmtree(directory_name)
    for i in xrange(3):
        root.after(1,wrapped_partial(root_closer, root))
    for i in xrange(20):
        root.quit()

if __name__ == "__main__":
    import shutil, tempfile
    directory_name = tempfile.mkdtemp()
    shutil.copy(resource_path("docs/corvids_readme.html"), directory_name)
    shutil.copy(resource_path("docs/Corvid.png"), directory_name)
    shutil.copy(resource_path("docs/sample_fig.png"), directory_name)

    root = Tk()
    root.protocol("WM_DELETE_WINDOW", wrapped_partial(root_closer, root, directory_name))
    gui = CoreGUI(root, directory_name)
    root.mainloop()
