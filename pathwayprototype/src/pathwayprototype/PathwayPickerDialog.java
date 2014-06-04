/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ruthgrace
 */

package pathwayprototype;

import java.net.URI;
import java.awt.Desktop;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
public class PathwayPickerDialog {
    private String[] pathwayNames;
    private String[] pathwayLinks;
    private JDialog dialog;
    private JPanel panel;
    private JComboBox pathwayList;
    private JButton okButton;
    private JLabel pathwayLabel;
    private JTextArea description;
    public PathwayPickerDialog(String[] pathwayNames, String[] htmlLinks) {
        
        this.pathwayNames = pathwayNames;
        this.pathwayLinks = htmlLinks;
        JDialog emptyParent = null;
        dialog = new JDialog(emptyParent,"Choose a pathway to display in your browser");
        panel = new JPanel(); // you can give it a layout
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        description = new JTextArea("The following pathways were found to contain genes containing the mutations given. Pathways are displayed as a network. Mutated genes are colored blue, and all other pathway components are colored black.");
        description.setEditable(false);
        panel.add(description);
        
        pathwayLabel = new JLabel("Please select a pathway to be displayed in your browser.");
        panel.add(pathwayLabel);
        
        final JComboBox pathwayList = new JComboBox(pathwayNames);
        panel.add(pathwayList);
        
        okButton = new JButton("OPEN");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.out.println("LINK CLICKED: "+pathwayLinks[pathwayList.getSelectedIndex()]+ " at index; "+pathwayList.getSelectedIndex());
                Path currentRelativePath = Paths.get("");
                String pathToProgram = currentRelativePath.toAbsolutePath().toString();
                String pathStart = "file://";
                String slash = "/";
                openLink(pathStart+pathToProgram+slash+pathwayLinks[pathwayList.getSelectedIndex()]);
                System.out.println("link open attempted");
            }
        
        });
        panel.add(okButton);
        
        dialog.getContentPane().add(panel);
        dialog.pack();
        dialog.setVisible(true);
        
    }
    private static void openLink(String htmlFileLink) {
        try {
            final URI uri = new URI(htmlFileLink);
            if (Desktop.isDesktopSupported()) {
              try {
                Desktop.getDesktop().browse(uri);
              } catch (IOException e) { /* TODO: error handling */ }
            } else { /* TODO: error handling */ }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
