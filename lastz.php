<html>
<head>
<style>
	.formdiv {
		background-color: #cccccd;
		padding: 1%;
		border-style: dashed;
		border-color: darkgray;
	}
</style>
</head>
<body>
	<div class="formdiv">
	<form id="lastzform" method="post" action="lastz.php">
		<div>
			<label for="lastzform">Lastz View Form:</label>
		</div>
		<label><b>Chromosome Options</b></label>
		<div>
			&nbsp;&nbsp;<label for="chr">Chromosome: </label><select name="chr">
				<?php
					$chr=array();
					for ($i=1; $i<23; $i++) {
						$chr[$i-1]=sprintf("chr%02d", $i);
					}
					$chr[22]="chrX";
					$chr[23]="chrY";
					for ($i=0; $i<24; $i++) {
						echo "<option value=\"", $chr[$i], "\"";
						if ($chr[$i] == $_POST['chr'])
							echo " selected";
						echo ">", $chr[$i], "</option>";
					}
	
				?>
			</select>&nbsp;
			<label for="start">Start: </label><input name="start" value="<?php echo $_POST['start']; ?>"/>&nbsp;
			<label for="end">End: </label><input name="end" value="<?php echo $_POST['end']; ?>"/>&nbsp;
			<label for="gaps">Hide Gaps: </label><input name="gaps" type="checkbox" value="1" 
				<?php if (isset($_POST['gaps']))
					echo "checked";
				?>
			/>&nbsp;
			<label for="gname">Hide Gaps Position: </label><input name="gname" type="checkbox" value="0" 
				<?php if (isset($_POST['gname']))
					echo "checked";
				?>
			/>&nbsp;
			<label for="pname">Not to Transform ID: </label><input name="pname" type="checkbox" value="0" 
				<?php if (isset($_POST['pname']))
					echo "checked";
				?>
			/>&nbsp;
			<label for="galgn">Aligned to Gaps: </label><input name="galgn" type="checkbox" value="1" 
				<?php if (isset($_POST['galgn']))
					echo "checked";
				?>
			/>&nbsp;
		</div>
		<label><b>Image Options: </b></label>
		<div>
			&nbsp;&nbsp;<label for="width">Width: </label><input name="width" value="<?php 
				if (isset($_POST['width'])) {
					echo $_POST['width']; 
				} else {
					echo 2048;
				}
			?>"/>&nbsp;
			<label for="filter">Format: </label><select name="filter">
				<?php
					$ft=array("png", "jpg", "gif", "tif", "gd");
					foreach ($ft as $f) {
						echo "<option value=\"", $f, "\"";
						if ($f == $_POST['filter'])
 							echo " selected";
                        echo ">", $f, "</option>";
					}
				?>
			</select>&nbsp;
		</div>
		<label><b>Other Options: </b></label>
		<div>
			&nbsp;&nbsp;<label for="step">Step: </label><input name="step" value="<?php echo $_POST['step']; ?>"/>&nbsp;
			<button type="submit" name="bward"><?php echo " << "; ?></button>&nbsp;
			<button type="submit" name="zoomin"><?php echo " + "; ?></button>&nbsp;
			<button type="submit" name="zoomout"><?php echo " - "; ?></button>&nbsp;
			<button type="submit" name="fward"><?php echo " >> "; ?></button>&nbsp;
		</div>
		<div align="right">
			<button type="submit">Create</button>&nbsp;
			<button type="reset">Reset</button>&nbsp;
		</div>
	</form>
	</div>
	<hr/>
	<img src="<?php 
		if (isset($_POST)) {
			$url="/cgi/lastz.pl?"; 
			foreach ($_POST as $k => $v) {
				$url.="$k=$v&";
			}
		} else {
			$url="";
		}
		echo $url;
	?>" type="image/png" width="100%" id="lastz"/>
</body>
</html>

