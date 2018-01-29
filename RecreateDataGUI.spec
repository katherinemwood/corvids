# -*- mode: python -*-

block_cipher = None


a = Analysis(['RecreateDataGUI.py'],
             pathex=['/Users/katherinewood/Documents/Lab/CORVIDS/CORVIDS_5'],
             binaries=[],
             datas=[('./docs/*.png', 'docs'),('./docs/corvids_readme.html', 'docs')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='RecreateDataGUI',
          debug=False,
          strip=False,
          upx=True,
          console=False , icon='Corvid.icns')
app = BUNDLE(exe,
             name='RecreateDataGUI.app',
             icon='./Corvid.icns',
             bundle_identifier=None,
             info_plist={
                'NSHighResolutionCapable': 'True'})
